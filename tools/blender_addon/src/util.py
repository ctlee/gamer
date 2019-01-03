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

import gamer.pygamer as g

# bpy.context.tool_settings.mesh_select_mode = (True, False, False)

# DEFINITIONS
UNSETID = 0
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
    if not markerLayer:
        markerLayer = obj.data.polygon_layers_int.new('marker')
    return markerLayer.data


def getSelectedMesh():
    """
    Returns the first selected mesh
    """

    obj = bpy.context.active_object
    if obj:
        if obj.type == 'MESH':
             # Auto select object if in EDIT mode
            # This prevents being 'selected' in EDIT mode but really unselected
            if obj.mode == 'EDIT':
                obj.select = True
            if obj.select:
                return (True, obj)
            else:
                return (False, "Object is not selected! Select a MESH to use this feature.")
        else:
            return (False, "Active object is not a MESH. Please select a MESH to use this feature.")
    else:
        return (False, "No active object! Please select a MESH to use this feature.")


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
    mesh.update()
    mesh.calc_normals()
    return mesh


def blenderToGamer(obj=None, check_for_vertex_selection=False, map_boundaries=False):
    """
    Convert object to GAMer mesh.

    check_for_vertex_selection: True if selected vertices should be checked
    map_boundaries: True if boundary values should be mapped to markers
                    instead of boundary_id
    """
    success = True
    # Get the selected mesh
    if not obj:
        success, obj = getSelectedMesh()
    if not success:
        # If error, forward string along
        return (success, obj)

    with ObjectMode():
        bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
        # Grab vertices
        vertices, selected_vertices = getMeshVertices(obj, get_selected_vertices = True)
        # Get world location and offset each vertex with this value
        # translation = obj.location # Not needed after transform.apply
        gmesh = g.SurfaceMesh()   # Init GAMer SurfaceMesh

        def addVertex(co, sel): # functor to addVertices
            gmesh.addVertex(
                co[0],   # x position
                co[1],   # y position
                co[2],   # z position
                0,                        # marker
                bool(sel)                 # selected flag
            )

        # Check we have vertices selected
        if check_for_vertex_selection and not selected_vertices:
            return (False, "No vertices are selected.")

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

        # Get list of faces
        faces = obj.data.polygons

        # Transfer data from blender mesh to gamer mesh
        for face in faces:
            vertices = deque(face.vertices)
            if len(vertices) != 3:
                return (False, "Encountered a non-triangular face. GAMer only works with triangulated meshes.")

            # Get the orientation from Blender
            max_val = max(vertices)
            max_idx = vertices.index(max_val)
            if max_idx != 2:
                vertices.rotate(2-max_idx)
            orientation = 1
            if(vertices[0] < vertices[1]):
                orientation = -1
            gmesh.insertFace(g.FaceKey(vertices), g.Face(orientation, boundaries[face.index], False))
    # Ensure all face orientations are set
    g.init_orientation(gmesh)
    g.check_orientation(gmesh)
    vol = g.getVolume(gmesh)
    print(vol)
    return (True, gmesh)


def gamerToBlender(gmesh,
        obj = None,
        mesh_name="gamer_improved"):
    # Check arguments
    if not isinstance(gmesh, g.SurfaceMesh):
        return (False, "gamerToBlender expected a SurfaceMesh")

    if not obj:
        success, obj = getSelectedMesh()
    if not success:
        # If error, forward string along
        return (success, obj)

    with ObjectMode():
        verts = []      # Array of vertex coordinates
        idxMap = {}     # Dictionary of gamer indices to renumbered indices
        un_selected_vertices = []
        for i, vid in enumerate(gmesh.vertexIDs()):
            v = vid.data()
            verts.append((v[0], v[1], v[2]))
            idxMap[gmesh.getName(vid)[0]] = i    # Reindex starting at 0
            if not v.selected:
                un_selected_vertices.append(i)

        faces = []
        markersList = []
        for i, fid in enumerate(gmesh.faceIDs()):
            fName = gmesh.getName(fid)
            face = fid.data()

            if face.orientation == 0:
                return (False, "gamerToBlender: Found face with undefined orientation")
            if face.orientation == -1:
                faces.append((idxMap[fName[0]], idxMap[fName[1]], idxMap[fName[2]]))
            else:
                faces.append((idxMap[fName[2]], idxMap[fName[1]], idxMap[fName[0]]))
            markersList.append(face.marker)

        # if create_new_mesh:
        #     # Create new object and mesh
        #     bmesh = createMesh(mesh_name, verts, faces)
        #     obj = bpy.data.objects.new(mesh_name, bmesh)

        #     bpy.context.scene.objects.link(obj) # link object to scene

        #     bpy.ops.object.select_all(action='DESELECT')
        #     obj.select=True

        #     ml = markers.getMarkerLayer(obj)
        #     for i, marker in enumerate(markersList):
        #         ml[i].value = marker
        #     markers.copyBoundaries(currObj, obj)
        #     currObj = obj
        #     print("Create_new_mesh = False")

        orig_mesh = obj.data
        bmesh = createMesh('gamer_tmp', verts, faces)
        obj.data = bmesh
        ml = getMarkerLayer(obj)
        for i, marker in enumerate(markersList):
            ml[i].value = marker

        bpy.data.meshes.remove(orig_mesh)
        bmesh.name = mesh_name


        def toggle(vert, val):
            vert = val

        [toggle(v.select,False) for v in bmesh.vertices if v.index in un_selected_vertices]
        obj.select = True

    # Repaint boundaries
    obj.gamer.repaint_boundaries(bpy.context)
    return (True, None)



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
        me = obj.to_mesh(depsgraph=bpy.context.depsgraph, apply_modifiers=True)
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

    # TODO. remove all customdata layers.
    # would save ram

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
