# ***************************************************************************
# This file is part of the GAMer software.
# Copyright (C) 2016-2018
# by Christopher Lee, Tom Bartol, John Moody, Rommie Amaro, J. Andrew McCammon,
#    and Michael Holst

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.

# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# ***************************************************************************

import bpy
import bmesh

import numpy as np
from collections import deque

import blendgamer.pygamer as pygamer
import blendgamer.pygamer.surfacemesh as sm

# bpy.context.tool_settings.mesh_select_mode = (True, False, False)

# DEFINITIONS
UNSETID = 0     # Default should match MeshPolygonIntProperty default
UNSETMARKER = -1
materialNamer = lambda bnd_id: "%s_mat"%(bnd_id)
##

class ObjectMode():
    def __enter__(self):
        self.mode = bpy.context.active_object.mode
        bpy.ops.object.mode_set(mode='OBJECT')

    def __exit__(self, type, value, traceback):
        # print("Changing to %s mode"%(self.mode))
        bpy.ops.object.mode_set(mode=self.mode)

def getBoundaryMaterial(boundary_id):
    mats = bpy.data.materials
    for mat in mats:
        if mat.gamer.boundary_id == boundary_id:
            return mat

def getMarkerLayer(obj):
    if obj.mode != 'OBJECT':
        raise RuntimeError("Blender Layers (Markers) can only be accessed in 'OBJECT' mode.")
    markerLayer = obj.data.polygon_layers_int.get('marker')
    # If the layer doesn't exist yet, create it!
    if not markerLayer:
        markerLayer = obj.data.polygon_layers_int.new('marker')
    return markerLayer.data


def getActiveMeshObject(report):
    """
    Returns the active object if it is a mesh or None
    """
    obj = bpy.context.active_object
    if obj:
        if obj.type == 'MESH':
            #  # Auto select object if in EDIT mode
            # # This prevents being 'selected' in EDIT mode but really unselected
            # if obj.mode == 'EDIT':
            #     obj.select = True
            # if obj.select:
            return obj
        else:
            report({'ERROR'}, "Active object is not a MESH. Please select a MESH to use this feature.")
    else:
        report({'ERROR'}, "No active object! Please select a MESH to use this feature.")
    return None


def getMeshVertices(obj, get_selected_vertices=False):
    """
    Get the vertices of mesh object
    """
    mesh = obj.data
    vertToVec = lambda v : [v[0], v[1], v[2]]

    vertices = [vertToVec(v.co) for v in mesh.vertices]

    if get_selected_vertices:
        selected_indices = [v.index for v in mesh.vertices if v.select and not v.hide]
        return vertices, selected_indices
    else:
        return vertices


def createMesh(mesh_name, verts, faces):
    """
    @brief      Creates a new blender mesh from arrays of vertex and face
    """
    mesh = bpy.data.meshes.new(mesh_name)
    mesh.from_pydata(verts, [], faces)
    # mesh.validate(verbose=True)
    mesh.update()
    mesh.calc_normals()
    return mesh


def blenderToGamer(report, obj=None, map_boundaries=False, autocorrect_normals=True):
    """
    Convert object to GAMer mesh.

    map_boundaries: True if boundary values should be mapped to markers
                    instead of boundary_id

    @Returns    gamer.SurfaceMesh object or None
    """
    # Get the selected mesh
    if not obj:
        obj = getActiveMeshObject(report)
        # Ensure object is good
        if not obj:
            return False

    with ObjectMode():
        # Grab vertices
        vertices, selected_vertices = getMeshVertices(obj, get_selected_vertices = True)
        # Get world location and offset each vertex with this value
        gmesh = sm.SurfaceMesh()   # Init GAMer SurfaceMesh

        def addVertex(co, sel): # functor to addVertices
            gmesh.addVertex(sm.Vertex(
                co[0],   # x position
                co[1],   # y position
                co[2],   # z position
                0,
                bool(sel)                 # selected flag
            ))

        # If all vertices are selected
        if len(selected_vertices) == len(vertices):
            selected_vertices = np.ones(len(vertices), dtype=bool)
        else:
            selection = np.zeros(len(vertices), dtype=bool)
            selection[selected_vertices] = 1    # Set selected to True
            selected_vertices = selection

        # Zip args and pass to addVertex functor
        [addVertex(*args) for args in zip(vertices, selected_vertices)]

        ml = getMarkerLayer(obj)
        # Transfer boundary information
        if map_boundaries:
            bdryMap = {UNSETID: UNSETMARKER}
            for bdry in obj.gamer.boundary_list:
                bdryMap[bdry.boundary_id] = bdry.marker
            boundaries = [bdryMap[item.value] for item in ml.values()]
        else:
            boundaries = [item.value for item in ml.values()]

        edges = obj.data.edges
        for edge in edges:
            gmesh.insertEdge(list(edge.vertices), sm.Edge(bool(edge.select)))

        # Get list of faces
        faces = obj.data.polygons

        # Transfer data from blender mesh to gamer mesh
        for face in faces:
            vertices = deque(face.vertices)
            if len(vertices) != 3:
                report({'ERROR'}, "Encountered a non-triangular face. GAMer only works with triangulated meshes.")
                return None

            # Get the orientation from Blender
            max_val = max(vertices)
            max_idx = vertices.index(max_val)
            if max_idx != 2:
                vertices.rotate(2-max_idx)
            orientation = -1    # 1 < 0 < 2
            if(vertices[0] < vertices[1]):
                # 0 < 1 < 2
                orientation = 1
            gmesh.insertFace(list(vertices), sm.Face(orientation, boundaries[face.index], bool(face.select)))
    # Ensure all face orientations are set
    gmesh.init_orientation()
    gmesh.check_orientation()
    vol = gmesh.getVolume()
    if vol < 0:
        if autocorrect_normals:
            gmesh.flipNormals()
        else:
            report({'ERROR'}, "Mesh has negative volume. Recompute normals to be outward facing.")
            return None
    return  gmesh


def gamerToBlender(report, gmesh,
        obj = None,
        mesh_name="gamer_improved"):
    # Check arguments
    if not isinstance(gmesh, sm.SurfaceMesh):
        report({'ERROR'}, "gamerToBlender expected a pygamer.surfacemesh.SurfaceMesh object")
        return False

    if not obj:
        obj = getActiveMeshObject(report)
        if not obj:
            return False

    mode = obj.mode

    verts = []      # Array of vertex coordinates
    idxMap = {}     # Dictionary of gamer indices to renumbered indices
    for i, vid in enumerate(gmesh.vertexIDs):
        v = vid.data()
        verts.append((v[0], v[1], v[2]))
        idxMap[gmesh.getName(vid)[0]] = i    # Reindex starting at 0

    faces = []
    selectedFaces = []
    markersList = []
    for i, fid in enumerate(gmesh.faceIDs):
        fName = gmesh.getName(fid)
        face = fid.data()

        if face.selected:
            selectedFaces.append(i)

        if face.orientation == 1:
            faces.append((idxMap[fName[0]], idxMap[fName[1]], idxMap[fName[2]]))
        else:
            faces.append((idxMap[fName[2]], idxMap[fName[1]], idxMap[fName[0]]))
        markersList.append(face.marker)

    with ObjectMode():
        # Save old mesh data...
        # Can't delete yet or object goes out of context
        oldmesh = obj.data
        # Create new mesh and delete old
        newmesh = createMesh(mesh_name, verts, faces)
        obj.data = newmesh
        bpy.data.meshes.remove(oldmesh)

        ml = getMarkerLayer(obj)
        for i, marker in enumerate(markersList):
            ml[i].value = marker

    # Repaint boundaries
    obj.gamer.repaint_boundaries(bpy.context)
    # Deselect all first
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')

    bm = bmesh_from_object(obj)
    bm.faces.ensure_lookup_table()
    for f in selectedFaces:
        bm.faces[f].select_set(True)
    bmesh_to_object(obj, bm)
    bpy.ops.object.mode_set(mode=mode)
    return True



## Following functions are from 3D Print Addon
def clean_float(text):
    # strip trailing zeros: 0.000 -> 0.0
    index = text.rfind(".")
    if index != -1:
        index += 2
        head, tail = text[:index], text[index:]
        tail = tail.rstrip("0")
        text = head + tail
    return text

def bmesh_copy_from_object(obj, transform=True, triangulate=True, apply_modifiers=False):
    """
    Returns a transformed, triangulated copy of the mesh
    """

    assert obj.type == 'MESH'

    if apply_modifiers and obj.modifiers:
        import bpy
        me = obj.to_mesh(bpy.context.scene, apply_modifiers=True, settings='PREVIEW')
        bm = bmesh.new()
        bm.from_mesh(me)
        bpy.data.meshes.remove(me)
        del bpy
    else:
        me = obj.data
        if obj.mode == 'EDIT':
            bm_orig = bmesh.from_edit_mesh(me)
            bm = bm_orig.copy()
        else:
            bm = bmesh.new()
            bm.from_mesh(me)

    if transform:
        bm.transform(obj.matrix_world)

    if triangulate:
        bmesh.ops.triangulate(bm, faces=bm.faces)

    return bm


def bmesh_from_object(obj):
    """
    Object/Edit Mode get mesh, use bmesh_to_object() to write back.
    """
    me = obj.data
    is_editmode = (obj.mode == 'EDIT')
    if is_editmode:
        bm = bmesh.from_edit_mesh(me)
    else:
        bm = bmesh.new()
        bm.from_mesh(me)
    return bm


def bmesh_to_object(obj, bm):
    """
    Object/Edit Mode update the object.
    """
    me = obj.data
    is_editmode = (obj.mode == 'EDIT')
    if is_editmode:
        bmesh.update_edit_mesh(me, True)
    else:
        bm.to_mesh(me)
    # grr... cause an update
    if me.vertices:
        me.vertices[0].co[0] = me.vertices[0].co[0]


def bmesh_calc_area(bm):
    """
    Calculate the surface area.
    """
    return sum(f.calc_area() for f in bm.faces)


def bmesh_check_self_intersect_object(obj):
    """
    Check if any faces self intersect

    returns an array of edge index values.
    """
    import array
    import mathutils

    if not obj.data.polygons:
        return array.array('i', ())

    bm = bmesh_copy_from_object(obj, transform=False, triangulate=False)
    tree = mathutils.bvhtree.BVHTree.FromBMesh(bm, epsilon=0.00001)
    overlap = tree.overlap(tree)
    faces_error = {i for i_pair in overlap for i in i_pair}

    return array.array('i', faces_error)
