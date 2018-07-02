import bpy
from bpy.props import BoolProperty, CollectionProperty, EnumProperty, \
    FloatProperty, FloatVectorProperty, IntProperty, IntVectorProperty, \
    PointerProperty, StringProperty, BoolVectorProperty
from bpy.app.handlers import persistent
import gamer.pygamer as g

from . import boundary_markers
from . import tetrahedralization

# python imports
import os
import numpy as np
import collections

# we use per module class registration/unregistration
def register():
    bpy.utils.register_module(__name__)


def unregister():
    bpy.utils.unregister_module(__name__)

@persistent
def gamer_load_post(context):
    """ Initialize GAMer add  """
    print ( "load post handler: gamer_load_post() called"
)
    if not context:
        context = bpy.context
    scn = bpy.context.scene
    if not scn.gamer.initialized:
      scn.gamer.init_properties()


def panel_select_callback (self, context):
    self.panel_select_callback(context)


class GAMER_OT_coarse_dense(bpy.types.Operator):
    bl_idname = "gamer.coarse_dense"
    bl_label = "Coarse Dense Tris"
    bl_description = "Decimate selected dense areas of the mesh"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.scene.gamer.mesh_improve_panel.coarse_dense(context)
        return {'FINISHED'}


class GAMER_OT_coarse_flat(bpy.types.Operator):
    bl_idname = "gamer.coarse_flat"
    bl_label = "Coarse Flat Tris"
    bl_description = "Decimate selected flat areas of the mesh"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.scene.gamer.mesh_improve_panel.coarse_flat(context)
        return {'FINISHED'}


class GAMER_OT_smooth(bpy.types.Operator):
    bl_idname = "gamer.smooth"
    bl_label = "Smooth Tris"
    bl_description = "Smooth selected vertices of the mesh"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.scene.gamer.mesh_improve_panel.smooth(context)
        return {'FINISHED'}


class GAMER_OT_normal_smooth(bpy.types.Operator):
    bl_idname = "gamer.normal_smooth"
    bl_label = "Normal Smooth Surf"
    bl_description = "Smooth facet normals of selected faces of the mesh"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.scene.gamer.mesh_improve_panel.normal_smooth(context)
        return {'FINISHED'}


class GAMerMeshImprovementPropertyGroup(bpy.types.PropertyGroup):
  dense_rate = FloatProperty(
      name="CD_Rate", default=2.5, min=0.001, max=4.0, precision=4,
      description="The rate for coarsening dense areas")
  dense_iter = IntProperty(
      name="CD_Iter", default=1, min=1, max=15,
      description="The number of iterations for coarsening dense areas")
  flat_rate = FloatProperty(
      name="CF_Rate", default=0.016, min=0.00001, max=0.5, precision=4,
      description="The rate for coarsening flat areas")
  max_min_angle = IntProperty(
      name="Max_Min_Angle", default=15, min=10, max=20,
      description="The maximal minumum angle for smoothing")
  smooth_iter = IntProperty(
      name="S_Iter", default=6, min=1, max=50,
      description="The number of iterations for coarsening dense areas")
  preserve_ridges = BoolProperty( name="Preserve ridges", default=False)
  new_mesh = BoolProperty( name="Create new mesh", default=False)

  def coarse_dense (self, context):
      print("Calling coarse_dense")
      gmesh, boundaries = blenderToGamer()
      gmesh.coarse_dense(rate=self.dense_rate, numiter=self.dense_iter)
      gamer_to_blender(gmesh, boundaries, create_new_mesh=self.new_mesh)

  def coarse_flat (self, context):
      print("Calling coarse_flat")
      gmesh, boundaries = blenderToGamer()
      gmesh.coarse_flat(rate=self.flat_rate)
      gamer_to_blender(gmesh, boundaries, create_new_mesh=self.new_mesh)

  def smooth (self, context):
      print("Calling smooth")
      gmesh, boundaries = blenderToGamer()
      g.smoothMesh(gmesh, max_min_angle=self.max_min_angle, max_iter=self.smooth_iter, preserve_ridges=self.preserve_ridges)
      # gamer_to_blender(gmesh, boundaries, create_new_mesh=self.new_mesh)

  def normal_smooth (self, context):
      print("Calling smooth")
      gmesh, boundaries = blenderToGamer()
      gmesh.normal_smooth()
      gamer_to_blender(gmesh, boundaries, create_new_mesh=self.new_mesh)

  def draw_layout (self, context, layout ):
      row = layout.row()
      col = row.column()
      col.operator("gamer.coarse_dense",icon="OUTLINER_OB_LATTICE")
      col = row.column()
      col.prop(self, "dense_rate" )
      col = row.column()
      col.prop(self, "dense_iter" )

      row = layout.row()
      col = row.column()
      col.operator("gamer.coarse_flat",icon="MOD_TRIANGULATE")
      col = row.column()
      col.prop(self, "flat_rate" )

      row = layout.row()
      col = row.column()
      col.operator("gamer.smooth",icon="OUTLINER_OB_MESH")
      col = row.column()
      col.prop(self, "max_min_angle" )
      col = row.column()
      col.prop(self, "smooth_iter" )

      row = layout.row()
      row.prop(self, "preserve_ridges" )

      row = layout.row()
      col = row.column()
      col.operator("gamer.normal_smooth",icon="SMOOTHCURVE")

#      row = layout.row()
#      row.prop(self, "new_mesh" )


  def draw_panel ( self, context, panel ):
      layout = panel.layout
      self.draw_layout ( context, layout )



class GAMerMainPanelPropertyGroup(bpy.types.PropertyGroup):
    mesh_improve_select = BoolProperty ( name="mesh_improve_sel", description="Surface Mesh Improvement", default=False, subtype='NONE', update=panel_select_callback)
    boundary_markers_select = BoolProperty ( name="boundary_markers_sel", description="Boundary Markers", default=False, subtype='NONE', update=panel_select_callback)
    tet_select = BoolProperty ( name="tet_sel", description="Tetrahedralization", default=False, subtype='NONE', update=panel_select_callback)
    select_multiple = BoolProperty ( name="select_multiple", description="Show Multiple Panels", default=False, subtype='NONE', update=panel_select_callback)
    last_state = BoolVectorProperty ( size=22 ) # Keeps track of previous button state to detect transitions

    def panel_select_callback ( self, context ):
        """
        Desired Logic:
          pin_state 0->1 with no others selected:
            Show All
          pin_state 0->1 with just 1 selected:
            No Change (continue showing the currently selected, and allow more)
          pin_state 0->1 with more than 1 selected ... should NOT happen because only one panel should show when pin_state is 0
            Illegal state
          pin_state 1->0 :
            Hide all panels ... always
        """

        prop_keys = [ 'mesh_improve_select', 'boundary_markers_select', 'tet_select', 'select_multiple' ]

        pin_state = False

        if self.get('select_multiple'):
            pin_state = (self['select_multiple'] != 0)
        old_pin_state = (self.last_state[prop_keys.index('select_multiple')] != 0)

        if (old_pin_state and (not pin_state)):
            # Pin has been removed, so hide all panels ... always
            # print ("Hiding all")
            for k in prop_keys:
                self.last_state[prop_keys.index(k)] = False
                self[k] = 0
            self.last_state[prop_keys.index('select_multiple')] = False

        elif ((not old_pin_state) and pin_state):
            # Pin has been pushed
            # Find out how many panels are currently shown
            num_panels_shown = 0
            for k in prop_keys:
                if k != 'select_multiple':
                    if self.get(k):
                        if self[k] != 0:
                            num_panels_shown += 1
            # Check for case where no panels are showing
            if num_panels_shown == 0:
                # print ("Showing all")
                # Show all panels
                for k in prop_keys:
                    if self.get(k):
                        self[k] = 1
                        self.last_state[prop_keys.index(k)] = False

            self.last_state[prop_keys.index('select_multiple')] = True

        else:
            # Pin state has not changed, so assume some other button has been toggled

            # Go through and find which one has changed to positive, setting all others to 0 if not pin_state
            for k in prop_keys:
                if self.get(k):
                    if (self[k] != 0) and (self.last_state[prop_keys.index(k)] == False):
                        self.last_state[prop_keys.index(k)] = True
                    else:
                        if not pin_state:
                            self.last_state[prop_keys.index(k)] = False
                            self[k] = 0


    def draw_self(self, context, layout):

        # Draw all the panel selection buttons with labels in 2 columns:

        brow = layout.row()
        bcol = brow.column()
        bcol.prop ( self, "mesh_improve_select", icon='MESH_ICOSPHERE', text="Surface Mesh Improvement" )
        bcol = brow.column()
        bcol.prop ( self, "boundary_markers_select", icon='TPAINT_HLT', text="Boundary Marking" )

        brow = layout.row()
        bcol = brow.column()
        bcol.prop ( self, "tet_select", icon='MOD_SKIN', text="Tetrahedralization" )
        bcol = brow.column()
        if self.select_multiple:
            bcol.prop ( self, "select_multiple", icon='PINNED', text="Show All / Multiple" )
        else:
            bcol.prop ( self, "select_multiple", icon='UNPINNED', text="Show All / Multiple" )

        # Draw each panel only if it is selected

        if self.mesh_improve_select:
            layout.box() # Use as a separator
            layout.label ( "Surface Mesh Improvement", icon='MESH_ICOSPHERE' )
            context.scene.gamer.mesh_improve_panel.draw_layout ( context, layout )

        if self.boundary_markers_select:
            layout.box() # Use as a separator
            layout.label ( "Boundary Marking", icon='TPAINT_HLT' )
            active_obj = context.active_object
            if active_obj:
              active_obj.gamer.draw_layout ( context, layout )

        if self.tet_select:
            layout.box() # Use as a separator
            layout.label ( "Tetrahedralization", icon='MOD_SKIN' )
            context.scene.gamer.tet_group.draw_layout ( context, layout )




class GAMER_PT_main_panel(bpy.types.Panel):
    bl_label = "GAMer: Geometry-preserving Adaptive Mesher"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    bl_category = "GAMer"

    @classmethod
    def poll(cls, context):
        return (context.scene is not None)

    '''
    def draw_header(self, context):
        # LOOK HERE!! This is where the icon is actually included in the panel layout!
        # The icon() method takes the image data-block in and returns an integer that
        # gets passed to the 'icon_value' argument of your label/prop constructor or
        # within a UIList subclass
        img = bpy.data.images.get('cellblender_icon')
        #could load multiple images and animate the icon too.
        #icons = [img for img in bpy.data.images if hasattr(img, "icon")]
        if img is not None:
            icon = self.layout.icon(img)
            self.layout.label(text="", icon_value=icon)
    '''

    def draw(self, context):
        context.scene.gamer.main_panel.draw_self(context,self.layout)


class GAMerPropertyGroup(bpy.types.PropertyGroup):
  initialized = BoolProperty(name="GAMer Initialized", default=False)
  gamer_version = StringProperty(name="GAMer Version", default="0")
  boundary_id_counter = IntProperty(name="GAMer Boundary id Counter")

  main_panel = PointerProperty(
    type=GAMerMainPanelPropertyGroup,
    name="GAMer Main Panel")

  mesh_improve_panel = PointerProperty(
    type=GAMerMeshImprovementPropertyGroup,
    name="GAMer Surface Mesh Improvement")

  tet_group = PointerProperty(
    type=tetrahedralization.GAMerTetrahedralizationPropertyGroup,
    name="GAMer Tetrahedralization")

  def allocate_boundary_id ( self ):
    self.boundary_id_counter += 1
    boundary_id = "bnd_id_%d" % (self.boundary_id_counter)
    return self.boundary_id_counter, boundary_id

  def init_properties ( self ):
    self.gamer_version = "0.1"
    self.boundary_id_counter = 0

    if not bpy.data.materials.get('bnd_unset_mat') :
      bnd_unset_mat = bpy.data.materials.new('bnd_unset_mat')
      bnd_unset_mat.use_fake_user = True
      bnd_unset_mat.gamer.boundary_id = 'bnd_unset'

    self.initialized = True



def setObjectMode(obj):
    editmode = obj.mode
    bpy.ops.object.mode_set(mode='OBJECT')
    return editmode


def setEditMode(obj):
    editmode = obj.mode
    bpy.ops.object.mode_set(mode='EDIT')
    return editmode


def restoreInteractionMode(obj,editmode):
    bpy.ops.object.mode_set(mode=editmode)


def vertToVec(v):
    return[v[0],v[1],v[2]]


def toggle(var,val):
    var = val


def getSelectedMesh(errorreport=True):
    "Returns the selected mesh"
    objs = bpy.context.selected_objects
    obj = objs[0] if len(objs) else None
    # myprint(obj)
    if not obj or obj.type != 'MESH':
        if errorreport:
            print("expected a selected mesh")
        return None
    return obj



def getMeshVertices(obj, selected=False):
    """
    @brief      { item_description }

    @param      obj       The object
    @param      selected  If selected vertices should be grabbed separately

    @return     Get Vertices
    """
    mesh = obj.data
    if selected:
        vert_indices = [v.index for v in mesh.vertices if v.select and not v.hide]
        vertices = [vertToVec(mesh.vertices[vi].co) for vi in vert_indices]
        return vertices, vert_indices
    else:
        vertices = [vertToVec(v.co) for v in mesh.vertices]
        return vertices


def getMeshFaces( obj, selected = False):
    mesh = obj.data
    if selected :
        mfaces_indices = [face.index for face in mesh.polygons
                         if face.select and not face.hide]
        faces = [mesh.polygons[fi].vertices for fi in mfaces_indices]
        return faces, mfaces_indices
    else :
        faces = [f.vertices for f in mesh.polygons]
        return faces


def getBoundaryFaces(boundary):
    if not "faces" in boundary:
        return []
    all_faces = []
    for faces in list(boundary["faces"].values()):
        all_faces.extend(faces)
    return all_faces


def setBoundaryFaces(boundary, faces):
    "Set faces in boundary props"
    if not "faces" in boundary:
        return
    # Maximal indices in a array prop in Blender is 32767
    max_ind = 32767
    num_sub_arrays = int(len(faces)/max_ind)+1

    # If the faces already exist delete it and re attach it
    if "faces" in boundary:
        for key in boundary["faces"]:
            del boundary["faces"][key]
        del boundary["faces"]

    boundary["faces"] = {}
    for ind in range(num_sub_arrays):
        boundary["faces"]["F%d"%ind] = faces[ind*max_ind: \
                                             min((ind+1)*max_ind, len(faces))]


def createMesh(mesh_name, verts, faces):
    """
    @brief      Creates a new blender mesh from arrays of vertex and f

    @param      mesh_name  The mesh name
    @param      verts      The vertexes
    @param      faces      The faces

    @return     New blender mesh object
    """
    mesh = bpy.data.meshes.new(mesh_name)
    mesh.from_pydata(verts, [], faces)
    mesh.update()
    mesh.calc_normals()
    return mesh



def blenderToGamer(obj=None, check_for_vertex_selection=True):
    """
    @brief      Transfer active mesh to GAMer format

    @param      obj                         The object
    @param      check_for_vertex_selection  The check for vertex selection

    @return     Gamer Mesh & Boundaries
    """
    # Get the selected mesh
    if obj is None:
        obj = getSelectedMesh()
    if obj is None:
        return None, None

    editmode = setObjectMode(obj) # Ensure edit mode is off
    # Grab vertices
    vertices, selected_vertices = getMeshVertices(obj, selected=True)
    # Get world location and offset each vertex with this value
    translation = obj.location
    gmesh = g.SurfaceMesh()   # Init Gamer SurfaceMesh

    def addVertex(co, sel): # functor to addVertices
        gmesh.addVertex(g.Vertex(co[0] + translation[0],
                     co[1] + translation[1],
                     co[2] + translation[2],
            0, bool(sel)))

    # Check we have vertices selected
    if check_for_vertex_selection and not selected_vertices:
        print("blenderToGamer: There are no selected vertices.")
        return None, None

    # If all vertices are selected
    if len(selected_vertices) == len(vertices):
        selected_vertices = np.ones(len(vertices), dtype=bool)
    else:
        selection = np.zeros(len(vertices), dtype=bool)
        selection[selected_vertices] = 1    # Set selected to True
        selected_vertices = selection
    # Zip args and pass to addVertex functor
    [addVertex(*args) for args in zip(vertices, selected_vertices)]

    # Transfer boundary information
    boundaries = obj.get('boundaries')
    if not boundaries:
        obj['boundaries'] = {}
        boundaries = obj['boundaries']  # TODO (9): This line might be redundant CTL

    # map of indices to marker
    indexMap = g.index_map(len(obj.data.polygons))
    fillmap = lambda boundary: indexMap.updatemap(boundary['marker'], boundary['faces'].values()[0].to_list())
    # Iterate over the faces and transfer marker information
    map(fillmap, boundaries.values())

    # Get list of faces
    faces = obj.data.polygons

    # Transfer data from blender mesh to gamer mesh
    for face in faces:
        vertices = collections.deque(face.vertices)
        if len(vertices) != 3:
            # self.drawError(errormsg="expected mesh with only triangles in")
            print("Error: encountered a non-triangular face. Expected a triangulated mesh.")
            # self.waitingCursor(0)
            restoreInteractionMode(obj,editmode)
            return None, None

        # Get the orientation from Blender
        max_val = max(vertices)
        max_idx = vertices.index(max_val)
        if max_idx != 2:
            vertices.rotate(2-max_idx)
        orientation = 1
        if(vertices[0] < vertices[1]):
            orientation = -1
        gmesh.insertFace(g.FaceKey(vertices), g.Face(orientation, indexMap[face.index], False))

    # Ensure all face orientations are set
    g.compute_orientation(gmesh)

    # Restore editmode
    restoreInteractionMode(obj,editmode)
    return gmesh, boundaries



def gamerToBlender(gmesh, create_new_mesh=False, mesh_name="gamer_improved"):
    #myprint("gamer_to_blender ",gmesh)
    # Check arguments
    if not isinstance(gmesh, g.SurfaceMesh):
    # self.drawError(errormsg="expected a GAMer SurfaceMesh")
        print("expected a SurfaceMesh")

    # Get scene
    scn = bpy.context.scene
    # self.waitingCursor(1)

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
    bdryDict = dict()   # Dictionary of marker to list of face indices
    for i, fid in enumerate(gmesh.faceIDs()):
        fName = gmesh.getName(fid)
        face = fid.data()

        if face.orientation == -1:
            faces.append((idxMap[fName[0]], idxMap[fName[1]], idxMap[fName[2]]))
        else:
            faces.append((idxMap[fName[2]], idxMap[fName[1]], idxMap[fName[0]]))
        if face.marker != -1:
            if face.marker not in bndryDict:
                bdryDict[face.marker] = [face.marker]
            else:
                bdryDict[face.marker].append(face.marker)

    if create_new_mesh:
        # Create new object and mesh
        bmesh = createMesh(mesh_name, verts, faces)
        obj = bpy.data.objects.new(mesh_name, bmesh)

        scn.objects.link(obj)
        bpy.ops.object.select_all(action='DESELECT')
        obj.select=True



        if boundaries:
          for bnd_id in boundaries.keys():
              boundary = boundaries[bnd_id]
              new_boundaries[bnd_id] = dict(marker=, faces={})

        # Set the property dictionary
        # FIXME: Is this safe? Is boundaries always something I can use?
        obj['boundaries'] = boundaries
    else:
        orig_mesh = obj.data
        bmesh = createMesh('gamer_tmp', verts, faces)
        obj.data = bmesh
        bpy.data.meshes.remove(orig_mesh)
        bmesh.name = mesh_name

    #myprint("un_selected_vertices ", len(un_selected_vertices))
    [toggle(v.select,False) for v in bmesh.vertices if v.index in un_selected_vertices]

    # Restore editmode
    restoreInteractionMode(obj,editmode)

    # Repaint boundaries
    obj.gamer.repaint_boundaries(bpy.context)

    # self.waitingCursor(0)
    # self.updateViewer()

