# This branch was created to work on the tetrahedralization.

import bpy
from bpy.props import BoolProperty, CollectionProperty, EnumProperty, \
    FloatProperty, FloatVectorProperty, IntProperty, IntVectorProperty, \
    PointerProperty, StringProperty, BoolVectorProperty
import mathutils
import gamer
import gamer_addon.gamer_gui

# python imports
import os
import numpy as np


# we use per module class registration/unregistration
def register():
    bpy.utils.register_module(__name__)


def unregister():
    bpy.utils.unregister_module(__name__)


# Tetrahedralization Operators:

class GAMER_OT_set_tet_path(bpy.types.Operator):
    bl_idname = "gamer.set_tet_path"
    bl_label = "Set Tetrahedralization Path"
    bl_description = ("Set Tetrahedralization Path")
    bl_options = {'REGISTER'}

    filepath = bpy.props.StringProperty(subtype='FILE_PATH', default="")

    def execute(self, context):
        context.scene.gamer.tet_group.tet_path = self.filepath
        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


class GAMER_OT_tet_domain_add(bpy.types.Operator):
    bl_idname = "gamer.tet_domain_add"
    bl_label = "Add a Tet Domain"
    bl_description = "Add a new tetrahedralization domain"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.scene.gamer.tet_group.add_tet_domain(context)
        return {'FINISHED'}

class GAMER_OT_tet_domain_remove(bpy.types.Operator):
    bl_idname = "gamer.tet_domain_remove"
    bl_label = "Remove a Tet Domain"
    bl_description = "Remove a tetrahedralization domain"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.scene.gamer.tet_group.remove_active_tet_domain(context)
        self.report({'INFO'}, "Deleted Active Tet Group")
        return {'FINISHED'}

class GAMER_OT_tet_domain_remove_all(bpy.types.Operator):
    bl_idname = "gamer.tet_domain_remove_all"
    bl_label = "Remove a Tet Domain"
    bl_description = "Remove a tetrahedralization domain"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.scene.gamer.tet_group.remove_all_tet_domains(context)
        self.report({'INFO'}, "Deleted Active Tet Group")
        return {'FINISHED'}


class GAMER_OT_generic_button(bpy.types.Operator):
    bl_idname = "gamer.generic_button"
    bl_label = "Generic Button"
    bl_description = ("Generic Button")
    bl_options = {'REGISTER'}

    def execute(self, context):
        print ( "Executed" )
        return {'FINISHED'}

    def invoke(self, context, event):
        return self.execute(context)


class GAMER_OT_tetrahedralize(bpy.types.Operator):
    bl_idname = "gamer.tetrahedralize"
    bl_label = "Tetrahedralize"
    bl_description = ("Tetrahedralize")
    bl_options = {'REGISTER'}

    def execute(self, context):
        context.scene.gamer.tet_group.tetrahedralize()
        return {'FINISHED'}

    def invoke(self, context, event):
        return self.execute(context)


class GAMerTetDomainPropertyGroup(bpy.types.PropertyGroup):
    # name = StringProperty()  # This is a reminder that "name" is already defined for all subclasses of PropertyGroup
    domain_id = IntProperty ( name="ID", default=-1, description="Domain ID" )
    object_name = StringProperty ( name="ObjName", default="", description="Object Name" )
    marker = IntProperty ( name="Marker", default=-1, description="Domain Marker Integer" )
    is_hole = BoolProperty ( name="Hole", default=False, description="Use this domain as a hole" )
    constrain_vol  = BoolProperty ( name="Constrain Volume", default=False, description="Constrain Volume" )
    vol_constraint = FloatProperty ( name="Vol Constraint", default=10.0, description="Volume Constraint" )
    
    def draw_layout ( self, layout ):
        row = layout.row()
        col = row.column()
        col.prop ( self, "is_hole", text="Use Domain as a Hole" )
        if not self.is_hole:
            col = row.column()
            col.prop ( self, "marker" )
            row = layout.row()
            col = row.column()
            col.prop ( self, "constrain_vol" )
            if self.constrain_vol:
                col = row.column()
                col.prop ( self, "vol_constraint" )

    def draw_item_in_row ( self, row ):
        col = row.column()
        col.label ( str(self.object_name) )
        col = row.column()
        col.label ( "Domain ID: " + str(self.domain_id) )
        col = row.column()
        if self.is_hole:
            col.label ( "Hole" )
        else:
            col.label ( "Domain Marker: " + str(self.marker) )


class GAMer_UL_domain(bpy.types.UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        # The draw_item function is called for each item of the collection that is visible in the list.
        #   data is the RNA object containing the collection,
        #   item is the current drawn item of the collection,
        #   icon is the "computed" icon for the item (as an integer, because some objects like materials or textures
        #        have custom icons ID, which are not available as enum items).
        #   active_data is the RNA object containing the active property for the collection (i.e. integer pointing to the
        #        active item of the collection).
        #   active_propname is the name of the active property (use 'getattr(active_data, active_propname)').
        #   index is index of the current item in the collection.
        #   flt_flag is the result of the filtering process for this item.
        #   Note: as index and flt_flag are optional arguments, you do not have to use/declare them here if you don't
        #         need them.

        # The item will be a GAMerTetDomainPropertyGroup
        # Let it draw itself in a new row:

        item.draw_item_in_row ( layout.row() )



# Callbacks for all Property updates appear to require global (non-member) functions.
# This is circumvented by simply calling the associated member function passed as self:

def check_formats_callback(self, context):
    self.check_formats_callback(context)
    return

class GAMerTetrahedralizationPropertyGroup(bpy.types.PropertyGroup):

  tet_path = StringProperty ( name="tet_path", default="", description="Path and file prefix for tetrahedralization output files" )

  generic_float = FloatProperty( name="Generic Float", default=123.456, min=0.0, max=1000, precision=4, description="A Generic Float Value")
  generic_int = IntProperty( name="Generic Int", default=5, min=1, max=10, description="A Generic Int Value")
  generic_boolean = BoolProperty( name="Generic Bool", default=False, description="A Generic Boolean Value")

  domain_list = CollectionProperty(type=GAMerTetDomainPropertyGroup, name="Domain List")
  active_domain_index = IntProperty(name="Active Domain Index", default=0)
  next_id = IntProperty(name="Counter for Unique Domain IDs", default=1)  # Start ID's at 1 to confirm initialization

  show_settings = BoolProperty( name="Tetrahedralization Settings", default=False, description="Show more detailed settings")

  min_dihedral = FloatProperty ( name="Min Dihedral", default=10.0, description="Minimum Dihedral in Degrees" )
  max_aspect_ratio = FloatProperty ( name="Max Aspect Ratio", default=1.3, description="Maximum Aspect Ratio" )

  ho_mesh = BoolProperty ( name="Higher order mesh generation", default=False, description="Higher order mesh generation" )

  dolfin = BoolProperty ( name="DOLFIN", default=False, description="Generate DOLFIN output", update=check_formats_callback )
  diffpack = BoolProperty ( name="Diffpack", default=False, description="Generate Diffpack output", update=check_formats_callback )
  carp = BoolProperty ( name="Carp", default=False, description="Generate Carp output", update=check_formats_callback )
  fetk = BoolProperty ( name="FEtk", default=False, description="Generate FEtk output", update=check_formats_callback )
  
  status = StringProperty ( name="status", default="" )

  def check_formats_callback(self, context):
      if self.dolfin or self.diffpack or self.carp or self.fetk:
          self.status = ""

  def draw_layout ( self, context, layout ):

      row = layout.row()
      row.label ( "Domains" )

      row = layout.row()
      col = row.column()

      col.template_list("GAMer_UL_domain", "",
                        self, "domain_list",
                        self, "active_domain_index",
                        rows=2)

      col = row.column(align=True)
      col.operator("gamer.tet_domain_add", icon='ZOOMIN', text="")
      col.operator("gamer.tet_domain_remove", icon='ZOOMOUT', text="")
      col.operator("gamer.tet_domain_remove_all", icon='X', text="")

      if len(self.domain_list) > 0:
          domain = self.domain_list[self.active_domain_index]

          row = layout.row()
          row.label ( "Active Index = " + str ( self.active_domain_index ) + ", ID = " + str ( domain.domain_id ) )
          
          domain.draw_layout ( layout )

          box = layout.box()
          row = box.row(align=True)
          row.alignment = 'LEFT'
          if not self.show_settings:
              row.prop(self, "show_settings", icon='TRIA_RIGHT', emboss=False)
          else:
              row.prop(self, "show_settings", icon='TRIA_DOWN', emboss=False)

              row = box.row()
              row.operator("gamer.set_tet_path", text="Set Output File Prefix", icon='FILESEL')
              row = box.row()
              row.label ( self.tet_path )

              row = box.row()
              col = row.column()
              col.prop ( self, "min_dihedral" )
              col = row.column()
              col.prop ( self, "max_aspect_ratio" )

              row = box.row()
              row.prop ( self, "ho_mesh" )

              row = box.row()
              row.label ( "Output Formats:" )

              row = box.row()
              sbox = row.box()

              row = sbox.row()
              col = row.column()
              col.prop ( self, "dolfin" )
              col = row.column()
              col.prop ( self, "diffpack" )

              row = sbox.row()
              col = row.column()
              col.prop ( self, "carp" )
              col = row.column()
              col.prop ( self, "fetk" )

          row = layout.row()
          icon = 'PROP_OFF'
          if self.dolfin or self.diffpack or self.carp or self.fetk:
              icon = 'COLOR_RED'
          row.operator ( "gamer.tetrahedralize", text="Tetrahedralize", icon=icon )
          if len(self.status) > 0:
              row = layout.row()
              row.label ( self.status, icon="ERROR" )


  def add_tet_domain ( self, context):
      print("Adding a Tet Domain")
      """ Add a new tet domain to the list of tet domains for each selected object """
      #mcell = context.scene.mcell

      # From the list of selected objects, only add MESH objects.
      objs = [obj for obj in context.selected_objects if obj.type == 'MESH']
      if len(objs) > 0:
          for obj in objs:
              # Check by name to see if it's already listed
              current_domain_names = [ d.object_name for d in self.domain_list ]
              print ( "Current domains = " + str(current_domain_names) )
              if not (obj.name in current_domain_names):
                  new_id = self.allocate_available_id()  # Do this first to check for empty list before adding
                  obj.gamer.include = True
                  new_dom = self.domain_list.add()
                  new_dom.domain_id = new_id
                  new_dom.marker = new_id
                  new_dom.object_name = obj.name
                  self.active_domain_index = len(self.domain_list)-1

  def remove_active_tet_domain ( self, context):
      print("Removing active Tet Domain")
      """ Remove the active tet domain from the list of domains """
      self.domain_list.remove ( self.active_domain_index )
      self.active_domain_index -= 1
      if self.active_domain_index < 0:
          self.active_domain_index = 0

  def remove_all_tet_domains ( self, context):
      print("Removing All Tet Domains")
      """ Remove all tet domains from the list of domains """
      while len(self.domain_list) > 0:
          self.domain_list.remove ( 0 )
      self.active_domain_index = 0


  def allocate_available_id ( self ):
      """ Return a unique domain ID for a new domain """
      print ( "Next ID is " + str(self.next_id) )
      if len(self.domain_list) <= 0:
          # Reset the ID to 1 when there are no more molecules
          self.next_id = 1
      self.next_id += 1
      return ( self.next_id - 1 )

  def draw_panel ( self, context, panel ):
      layout = panel.layout
      self.draw_layout ( context, layout )

  def tetrahedralize ( self ):
      print ( "######################## Begin Tetrahedralize ########################" )

      filename = self.tet_path
      if not (self.dolfin or self.diffpack or self.carp or self.fetk):
          self.status = "Please select an output format in Tetrahedralization Settings"
          print ( self.status )
      else:
          self.status = ""
          mesh_formats = []

          if self.dolfin:
              mesh_formats.append("dolfin")
          if self.diffpack:
              mesh_formats.append("diffpack")
          if self.carp:
              mesh_formats.append("carp")
          if self.fetk:
              mesh_formats.append("mcsf")

          # What is this? It doesn't show up in the GAMer panel.
          #if self.getVal(self.tetparams["mcsf_format"]):
          #    mesh_formats.append("mcsf")

          # Get gamer mesh
          gmeshes = []
          boundary_markers = []

          for (obj_name,tet_domain) in [ (d.object_name,d) for d in self.domain_list ]:
              print ( "obj_name = " + obj_name + ", tet_domain = " + str(tet_domain) )

          current_domain_names = [ d.object_name for d in self.domain_list ]
          print ( "Current domains = " + str(current_domain_names) )
          for obj_name in current_domain_names:
              obj = bpy.data.objects[obj_name]
              if obj.gamer.include:
                  gmesh, boundaries = gamer_addon.gamer_gui.blender_to_gamer(obj=obj, create_new_mesh=False, check_for_vertex_selection=True)
                  if gmesh == None:
                      print ( "blender_to_gamer returned a gmesh of None" )
                  else:
                      # Collect boundary information
                      for boundary_name, boundary in zip(boundaries.keys(), boundaries.values()):
                          boundary_markers.append((boundary["marker"], boundary_name))

                      print ("\nMesh %s: num verts: %d numfaces: %d" % (obj_name, gmesh.num_vertices, gmesh.num_faces))

                      # Set the domain data on the SurfaceMesh these are the per/domain items as_hole, marker, and volume constraints
                      gmesh.as_hole = tet_domain.is_hole
                      gmesh.marker = tet_domain.marker
                      gmesh.use_volume_constraint = tet_domain.constrain_vol
                      gmesh.volume_constraint = tet_domain.vol_constraint

                      # Write surface mesh to file for debug
                      gmesh.write_off("surfmesh_%s.off" % obj_name)

                      # Add the mesh
                      gmeshes.append(gmesh)


          # Tetrahedralize mesh
          if len(gmeshes) > 0:

              quality_str = "q%.1fqq%.1faA"%(self.max_aspect_ratio,self.min_dihedral)

              quality_str += "o2" if self.ho_mesh else ""

              print("TetGen quality string: " + quality_str)

              # Do the tetrahedralization

              gem_mesh = gamer.GemMesh(gmeshes, quality_str)

              # Store mesh to files

              tetmesh_formats  = ["dolfin", "mcsf", "diffpack",  "carp"]
              tetmesh_suffices = [  ".xml",   ".m",    ".grid", ".carp"]

              for fmt in mesh_formats:

                  try:

                      suffix = tetmesh_suffices[tetmesh_formats.index(fmt)]
                      print ( "Writing to " + fmt + " file: " + filename + suffix )

                      # If the format is diffpack or carp we need to add boundary names
                      if fmt in ["diffpack", "carp"]:
                          boundary_names = [b[1] for b in sorted(boundary_markers)]
                          print ( "Boundary names = " + str(boundary_names) )
                          ### This doesn't work for some reason ... RuntimeError: Expected a list of strings as second argument.
                          getattr(gem_mesh, "write_%s"%fmt)(filename+suffix, boundary_names)
                      else:
                          getattr(gem_mesh, "write_%s"%fmt)(filename+suffix)

                  except Exception as ex:

                      print ( "Error: Unable to write to " + fmt + " file: " + filename + suffix )
                      print ( "   " + str(ex) )

      print ( "######################## End Tetrahedralize ########################" )



"""
    # Unmodified code from gamer/swig/src/upy_gui.py

    def tetrahedralize_action(self, filename):
        "Callback function for the tetrahedralize File chooser"

        # Load options and materials from registry
        self.load_from_registry()

        # Grab all selected mesh formats
        mesh_formats = []
        if self.getVal(self.tetparams["dolfin_format"]):
            mesh_formats.append("dolfin")
        if self.getVal(self.tetparams["diffpack_format"]):
            mesh_formats.append("diffpack")
        if self.getVal(self.tetparams["carp_format"]):
            mesh_formats.append("carp")
        if self.getVal(self.tetparams["mcsf_format"]):
            mesh_formats.append("mcsf")

        # Get gamer mesh
        gmeshes = []
        boundary_markers = []
        for i, (name, domain) in enumerate(self.domains.items()):
            obj = self.helper.getObject(name)
            if obj is None:
                self.drawError(errormsg="The domain: '%s' is not a mesh in "\
                               "this scene" % name)

            gmesh, boundaries = self.host_to_gamer(obj, False)
            if gmesh is None:
                return

            # Collect boundary information
            for boundary_name, boundary in zip(boundaries.keys(), boundaries.values()):
                boundary_markers.append((boundary["marker"], boundary_name))

            myprint("\nMesh %d: num verts: %d numfaces: %d" %
                    (i, gmesh.num_vertices, gmesh.num_faces))

            # Set the domain data on the SurfaceMesh
            for name, value in domain.items():
                setattr(gmesh, name, value)
                myprint("%s : %d" %(name, int(value)))

            # Write surface mesh to file for debug
            gmesh.write_off("surfmesh%d.off" % i)

            # Add the mesh
            gmeshes.append(gmesh)

        #obj = self._get_selected_mesh(False)
        #if obj is None:
        #    return

        # Tetrahedralize mesh
        quality_str = "q%.1fqq%.1faA"%(self.getVal(self.tetparams["aspect_ratio"]),
                                       self.getVal(self.tetparams["dihedral_angle"]))

        quality_str += "o2" if self.getVal(self.tetparams["higher_order"]) else ""

        myprint("TetGen quality string:", quality_str)

        # Do the tetrahedralization
        self.waitingCursor(1)
        gem_mesh = gamer.GemMesh(gmeshes, quality_str)

        # Store mesh to files
        for format_ in mesh_formats:

            suffix = self.tetmesh_suffices[self.tetmesh_formats.index(format_)]

            # If the format is diffpack or carp we need to add boundary names
            if format_ in ["diffpack", "carp"]:
                boundary_names = [b[1] for b in sorted(boundary_markers)]
                getattr(gem_mesh, "write_%s"%format_)(filename+suffix, boundary_names)
            else:
                getattr(gem_mesh, "write_%s"%format_)(filename+suffix)

        self.waitingCursor(0)

"""

