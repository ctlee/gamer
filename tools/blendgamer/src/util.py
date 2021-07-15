# ***************************************************************************
# This file is part of the GAMer software.
# Copyright (C) 2016-2021
# by Christopher T. Lee and contributors
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, see <http:#www.gnu.org/licenses/>
# or write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# **************************************************************************

import bpy
import bmesh

import numpy as np
from collections import deque
from contextlib import contextmanager

import blendgamer.pygamer as pygamer
import blendgamer.pygamer.surfacemesh as sm

# DEFINITIONS
UNSETID = 0  # ID value for unset markers
# Default should match MeshPolygonIntProperty default
UNSETMARKER = -1  # Marker value to output for unset faces
materialNamer = lambda bnd_id: "%s_mat" % (bnd_id)
##


class ObjectMode:
    """Class defining a Blender Object mode context for with statement semantics

    Attributes
    ----------
    mode : str
        Enum of Blender's mode
    """

    def __enter__(self):
        """
        Enter runtime context to cache current mode and switch to Object mode.
        """
        self.mode = bpy.context.active_object.mode
        bpy.ops.object.mode_set(mode="OBJECT")

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Exit runtime context and restore mode prior to enter.
        """
        # print("Changing to %s mode"%(self.mode))
        bpy.ops.object.mode_set(mode=self.mode)
        # Uncomment to handle exceptions
        # return True


@contextmanager
def BMeshContext(obj):
    if obj.type != "MESH":
        raise RuntimeError("Expected an object of MESH type.")

    me = obj.data
    if obj.mode == "EDIT":
        bm = bmesh.from_edit_mesh(me)
    else:
        bm = bmesh.new()
        bm.from_mesh(me)

    yield bm

    if obj.mode == "EDIT":
        if bpy.app.version < (2, 80, 0):
            bmesh.update_edit_mesh(me, tessface=True)
        else:
            bmesh.update_edit_mesh(me, loop_triangles=True)
        # BMesh from edit meshes should never be freed!
        # https://developer.blender.org/T39121
    else:
        bm.to_mesh(me)
        bm.free()
        del bm


@contextmanager
def copiedBMeshContext(obj, transform=False, triangulate=False, apply_modifiers=False):
    """Context for a bmesh copied from the mesh. Changes are not written back to the original mesh."""
    if obj.type != "MESH":
        raise RuntimeError("Expected object of MESH type.")

    if apply_modifiers and obj.modifiers:
        if bpy.app.version < (2, 80, 0):
            me = obj.to_mesh(
                bpy.context.scene, apply_modifiers=True, settings="PREVIEW"
            )
        else:
            me = obj.to_mesh()
        bm = bmesh.new()
        bm.from_mesh(me)
        bpy.data.meshes.remove(me)
    else:
        me = obj.data
        if obj.mode == "EDIT":
            bm_orig = bmesh.from_edit_mesh(me)
            bm = bm_orig.copy()
        else:
            bm = bmesh.new()
            bm.from_mesh(me)

    if transform:
        bm.transform(obj.matrix_world)

    if triangulate:
        bmesh.ops.triangulate(bm, faces=bm.faces)

    yield bm

    bm.free()
    del bm


def get_material_by_bnd_id(boundary_id):
    """
    Gets a boundary material from bpy.data.materials by boundary ID.

    Parameters
    ----------
    boundary_id : int
        The boundary ID to search for.

    Returns
    -------
    bpy.types.Material or None
        Returns the material if found or None if not.
    """
    mats = bpy.data.materials
    for mat in mats:
        if mat.gamer.boundary_id == boundary_id:
            return mat
    return None


def getBndUnsetMat():
    """
    Get the material for unset boundaries.

    Initializes the material if it does not already exist.

    Returns
    -------
    bpy.types.Material
        Returns the material corresponding to unmarked boundaries.
    """
    bnd_unset_mat = get_material_by_bnd_id(UNSETID)
    # if bnd_unset_mat is not defined, then create it
    if bnd_unset_mat is None:
        bnd_unset_mat = bpy.data.materials.new("bnd_unset_mat")
        bnd_unset_mat.use_fake_user = True
        bnd_unset_mat.gamer.boundary_id = UNSETID

    return bnd_unset_mat


def getMarkerLayer(obj):
    """Gets the boundary marker layer data of an object

    Parameters
    ----------
    obj : bpy.types.Object
        Object of interest

    Returns
    -------
    bpy.types.bpy_prop_collection
        Containing `bpy.types.MeshPolygonIntProperty` values corresponding to marker boundaryIDs

    Raises
    ------
    RuntimeError
        Data layers can only be accessed in 'OBJECT' mode
    """
    if obj.mode != "OBJECT":
        raise RuntimeError(
            "Blender Layers (Markers) can only be accessed in 'OBJECT' mode."
        )
    markerLayer = obj.data.polygon_layers_int.get("marker")
    # If the layer doesn't exist yet, create it!
    if not markerLayer:
        markerLayer = obj.data.polygon_layers_int.new(name="marker")
    return markerLayer.data


def getBMeshMarkerLayer(bm):
    """Get the boundary marker layer from a BMesh

    Args:
        bm (bpy.BMesh): BMesh object

    Returns:
        MarkerLayer: Marker layer
    """
    markerLayer = bm.faces.layers.int.get("marker")
    if not markerLayer:
        markerLayer = bm.faces.layers.int.new("marker")
    return markerLayer


def getCurvatureLayer(obj, algo, curvatureType):
    """Get the data layer corresponding to a specific curvature type

    Parameters
    ----------
    obj : bpy.types.Object
        Object of interest
    algo : TYPE
        Description
    curvatureType : TYPE
        Description

    Returns
    -------
    TYPE
        Description

    Raises
    ------
    RuntimeError
    Description
    """
    name = "%s%s" % (algo, curvatureType)
    if obj.mode != "OBJECT":
        raise RuntimeError(
            "Blender Layers (%s) can only be accessed in 'OBJECT' mode." % (name)
        )
    layer = obj.data.vertex_layers_float.get(name)
    # If the layer doesn't exist yet, create it!
    if not layer:
        layer = obj.data.vertex_layers_float.new(name=name)
    return layer.data


def getActiveMeshObject():
    """
    Returns the active object if it is a mesh or None

    Returns
    -------
    bpy.types.Object
        Mesh object of interest

    Raises
    ------
    RuntimeError
        If active object is not a Mesh or no active object exists.
    """
    obj = bpy.context.active_object
    if obj:
        if obj.type == "MESH":
            #  # Auto select object if in EDIT mode
            # # This prevents being 'selected' in EDIT mode but really unselected
            # if obj.mode == 'EDIT':
            #     obj.select = True
            # if obj.select:
            return obj
        else:
            raise RuntimeError(
                "Active object is not a MESH. Please select a MESH to use this feature."
            )
    else:
        raise RuntimeError(
            "No active object! Please select a MESH to use this feature."
        )


def getMeshVertices(obj, get_selected_vertices=False):
    """
    Get the vertices of mesh object

    Parameters
    ----------
    obj : TYPE
        Description
    get_selected_vertices : bool, optional
        Description

    Returns
    -------
    TYPE
        Description
    """
    mesh = obj.data
    vertToVec = lambda v: [v[0], v[1], v[2]]

    vertices = [vertToVec(v.co) for v in mesh.vertices]

    if get_selected_vertices:
        selected_indices = [v.index for v in mesh.vertices if v.select and not v.hide]
        return vertices, selected_indices
    else:
        return vertices


def blender_to_gamer(obj=None, map_boundaries=False, autocorrect_normals=True):
    """Construct a GAMer mesh from a triangulated blender mesh

    Args:
        obj (, optional): Object to convert. Defaults to None which gets the active mesh object.
        map_boundaries (bool, optional): True if boundary values should be mapped to markers instead of boundary_id. Defaults to False.
        autocorrect_normals (bool, optional): Automatically flip normals so that volume is positive. Defaults to True.

    Raises:
        RuntimeError: Complains if something prevents conversion to a valid GAMer mesh

    Returns:
        gamer.SurfaceMesh: GAMer surface mesh from Blender mesh
    """
    # Get the selected mesh
    if not obj:
        obj = getActiveMeshObject()
        # Ensure object is good
        if not obj:
            return False

    with ObjectMode():
        # Grab vertices
        vertices, selected_vertices = getMeshVertices(obj, get_selected_vertices=True)
        # Get world location and offset each vertex with this value
        gmesh = sm.SurfaceMesh()  # Init GAMer SurfaceMesh

        def addVertex(co, sel):  # functor to addVertices
            gmesh.addVertex(
                sm.Vertex(
                    co[0],  # x position
                    co[1],  # y position
                    co[2],  # z position
                    0,
                    bool(sel),  # selected flag
                )
            )

        # If all vertices are selected
        if len(selected_vertices) == len(vertices):
            selected_vertices = np.ones(len(vertices), dtype=bool)
        else:
            selection = np.zeros(len(vertices), dtype=bool)
            selection[selected_vertices] = 1  # Set selected to True
            selected_vertices = selection

        # Zip args and pass to addVertex functor
        [addVertex(*args) for args in zip(vertices, selected_vertices)]

        ml = getMarkerLayer(obj)
        # Transfer boundary information
        if map_boundaries:
            bdryMap = {UNSETID: UNSETMARKER}
            for bdry in obj.gamer.markers.boundary_list:
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
                raise RuntimeError(
                    "Encountered a non-triangular face. GAMer only works with triangulated meshes."
                )

            # Get the orientation from Blender
            max_val = max(vertices)
            max_idx = vertices.index(max_val)
            if max_idx != 2:
                vertices.rotate(2 - max_idx)
            orientation = -1  # 1 < 0 < 2
            if vertices[0] < vertices[1]:
                # 0 < 1 < 2
                orientation = 1
            gmesh.insertFace(
                list(vertices),
                sm.Face(orientation, boundaries[face.index], bool(face.select)),
            )
    # Ensure all face orientations are set
    gmesh.init_orientation()
    gmesh.check_orientation()
    vol = gmesh.getVolume()
    if vol < 0:
        if autocorrect_normals:
            gmesh.flipNormals()
        else:
            raise RuntimeError(
                "Mesh has negative volume. Recompute normals to be outward facing."
            )
    return gmesh


def gamer_to_blender(gmesh, obj=None, mesh_name="gamer_improved"):
    """Pass GAMer Surface Mesh back to Blender]

    Args:
        gmesh (gamer.SurfaceMesh): GAMer surface mesh
        obj (bpy.types.Object, optional):  Object to operate on. Active object is used if unset. Defaults to None.
        mesh_name (str, optional): Name of the new mesh data block. Defaults to "gamer_improved".

    Raises:
        RuntimeError: If an undefined behavior encountered.
    """

    # Check arguments
    if not isinstance(gmesh, sm.SurfaceMesh):
        raise RuntimeError(
            "gamer_to_blender expected a pygamer.surfacemesh.SurfaceMesh object."
        )

    if not obj:
        obj = getActiveMeshObject()
        if not obj:
            raise RuntimeError("Active object is not a MESH.")

    mode = obj.mode

    verts = []  # Array of vertex coordinates
    idxMap = {}  # Dictionary of gamer indices to renumbered indices
    for i, vid in enumerate(gmesh.vertexIDs):
        v = vid.data()
        verts.append((v[0], v[1], v[2]))
        idxMap[gmesh.getName(vid)[0]] = i  # Reindex starting at 0

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
        if bpy.app.version < (2, 81, 0):
            # Save old mesh data...
            # Can't delete yet or object goes out of context
            oldmesh = obj.data
            # Create new mesh and delete old
            mesh = bpy.data.meshes.new(mesh_name)
            mesh.from_pydata(verts, [], faces)
            # mesh.validate(verbose=True)
            mesh.update()
            # mesh.calc_normals()
            obj.data = mesh
            bpy.data.meshes.remove(oldmesh)
        else:
            mesh = obj.data
            mesh.clear_geometry()
            mesh.from_pydata(verts, [], faces)
            mesh.update()
            # mesh.calc_normals()

        ml = getMarkerLayer(obj)
        for i, marker in enumerate(markersList):
            ml[i].value = marker

    # Repaint boundaries
    obj.gamer.markers.repaint_boundaries(bpy.context)
    # Deselect all first
    bpy.ops.object.mode_set(mode="EDIT")
    bpy.ops.mesh.select_all(action="DESELECT")
    bpy.ops.object.mode_set(mode="OBJECT")

    bm = bmesh_from_object(obj)

    with BMeshContext(obj) as bm:
        bm.faces.ensure_lookup_table()
        for f in selectedFaces:
            bm.faces[f].select_set(True)

    bpy.ops.object.mode_set(mode=mode)


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


def bmesh_from_object(obj):
    """
    Object/Edit Mode get mesh, use bmesh_to_object() to write back.

    Parameters
    ----------
    obj : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    me = obj.data
    is_editmode = obj.mode == "EDIT"
    if is_editmode:
        bm = bmesh.from_edit_mesh(me)
    else:
        bm = bmesh.new()
        bm.from_mesh(me)
    return bm


def bmesh_to_object(obj, bm):
    """
    Object/Edit Mode update the object.
    NOTE that this free's the bmesh

    Parameters
    ----------
    obj : TYPE
        Description
    bm : TYPE
        Description
    """
    me = obj.data
    is_editmode = obj.mode == "EDIT"
    if is_editmode:
        bmesh.update_edit_mesh(me, True)
    else:
        bm.to_mesh(me)
        bm.free()
    # grr... cause an update
    if me.vertices:
        me.vertices[0].co[0] = me.vertices[0].co[0]


def make_annotations(cls):
    """Converts class fields to annotations if running with Blender 2.8

    This is a helper function allow code from 2.79 to work with 2.8x versions
    of Blender. All calls in register and unregister to pass through this.

    Returns
    -------
    class
        Converted class
    """
    if bpy.app.version < (2, 80):
        return cls

    if bpy.app.version >= (2, 93, 0):
        bl_props = {
            k: v
            for k, v in cls.__dict__.items()
            if isinstance(v, bpy.props._PropertyDeferred)
        }
    else:
        bl_props = {k: v for k, v in cls.__dict__.items() if isinstance(v, tuple)}

    if bl_props:
        if "__annotations__" not in cls.__dict__:
            setattr(cls, "__annotations__", {})
        annotations = cls.__dict__["__annotations__"]
        for k, v in bl_props.items():
            annotations[k] = v
            delattr(cls, k)
    return cls
