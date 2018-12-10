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
from bpy.props import BoolProperty, CollectionProperty, EnumProperty, \
    FloatProperty, FloatVectorProperty, IntProperty, IntVectorProperty, \
    PointerProperty, StringProperty, BoolVectorProperty
from bpy.app.handlers import persistent
import gamer.pygamer as g

from . import util
ObjectMode = util.ObjectMode

# python imports
import os
import re
import copy
from collections import deque

UNSETID = 0
UNSETMARKER = -1

# we use per module class registration/unregistration
def register():
    bpy.utils.register_module(__name__)


def unregister():
    bpy.utils.unregister_module(__name__)


@persistent
def boundary_markers_load_post(context):
    pass

# Object Boundary Marker Operators:
class GAMER_OT_add_boundary(bpy.types.Operator):
    bl_idname = "gamer.add_boundary"
    bl_label = "Add New Boundary"
    bl_description = "Add a new boundary to an object"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.object.gamer.add_boundary(context)
        return {'FINISHED'}


class GAMER_OT_remove_boundary(bpy.types.Operator):
    bl_idname = "gamer.remove_boundary"
    bl_label = "Remove Boundary"
    bl_description = "Remove selected boundary from object"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.object.gamer.remove_boundary(context)
        return {'FINISHED'}


class GAMER_OT_remove_all_boundaries(bpy.types.Operator):
    bl_idname = "gamer.remove_all_boundaries"
    bl_label = "Remove All Boundaries"
    bl_description = "Remove all boundaries from object"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.object.gamer.remove_all_boundaries(context)
        return {'FINISHED'}


class GAMER_OT_assign_boundary_faces(bpy.types.Operator):
    bl_idname = "gamer.assign_boundary_faces"
    bl_label = "Assign Selected Faces To Boundary"
    bl_description = "Assign selected faces to boundary"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        bnd = context.object.gamer.get_active_boundary()
        if bnd:
            bnd.assign_boundary_faces(context)
        return {'FINISHED'}


class GAMER_OT_remove_boundary_faces(bpy.types.Operator):
    bl_idname = "gamer.remove_boundary_faces"
    bl_label = "Remove Selected Faces From Boundary"
    bl_description = "Remove selected faces from boundary"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        bnd = context.object.gamer.get_active_boundary()
        if bnd:
            bnd.remove_boundary_faces(context)
        return {'FINISHED'}


class GAMER_OT_select_boundary_faces(bpy.types.Operator):
    bl_idname = "gamer.select_boundary_faces"
    bl_label = "Select Faces of Selected Boundary"
    bl_description = "Select faces of selected boundary"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        bnd = context.object.gamer.get_active_boundary()
        if bnd:
            bnd.select_boundary_faces(context)
        return {'FINISHED'}


class GAMER_OT_deselect_boundary_faces(bpy.types.Operator):
    bl_idname = "gamer.deselect_boundary_faces"
    bl_label = "Deselect Faces of Selected Boundary"
    bl_description = "Deselect faces of selected boundary"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        bnd = context.object.gamer.get_active_boundary()
        if bnd:
            bnd.deselect_boundary_faces(context)
        return {'FINISHED'}

class GAMER_OT_select_all_boundary_faces(bpy.types.Operator):
    bl_idname = "gamer.select_all_boundary_faces"
    bl_label = "Select All Faces of Selected Boundary"
    bl_description = "Select all faces of selected boundary"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        for bnd in context.object.gamer.boundary_list:
            bnd.select_boundary_faces(context)
        return {'FINISHED'}


# Object Boundary Panel:
class GAMER_UL_check_boundary(bpy.types.UIList):
    def draw_item(self, context, layout, data, item, icon, active_data,
                  active_propname, index):
        if item.status:
            layout.label(item.status, icon='ERROR')
        else:
            layout.label(item.name, icon='FILE_TICK')

        split = layout.split(percentage = 0.5)
        col = split.column()
        col = split.column()
        mats = bpy.data.materials
        bnd_mat = [ mat for mat in mats if mat.gamer.boundary_id == item.boundary_id ][0]
        col.prop(bnd_mat, 'diffuse_color', text='')


# Boundary Callbacks:
#
def boundary_name_update(self, context):
    # print('Boundary_name_update')
    # active_bnd_index = context.object.gamer.active_bnd_index
    # bnd = context.object.gamer.boundary_list[active_bnd_index]
    context.object.gamer.boundary_name_update(context)


def boundary_marker_update(self, context):
    active_bnd_index = context.object.gamer.active_bnd_index
    bnd = context.object.gamer.boundary_list[active_bnd_index]
    bnd.boundary_marker_update(context)


def getMarkerLayer(obj):
    if obj.mode != 'OBJECT':
        raise RuntimeError("Blender Layers (Markers) can only be accessed in 'OBJECT' mode.")
    markerLayer = obj.data.polygon_layers_int.get('marker')
    if not markerLayer:
        markerLayer = obj.data.polygon_layers_int.new('marker')
    return markerLayer.data


class GAMerBoundaryMaterial(bpy.types.PropertyGroup):
    """
    @brief      Class for GAMer boundary material property group.
    """
    boundary_id = IntProperty(name="Unique ID of This Boundary", default = -1)

class GAMerBoundaryMarker(bpy.types.PropertyGroup):
    """
    @brief      Class for GAMer boundary markers property group.
    """
    name = StringProperty(name="Boundary Name", default="Boundary", update=boundary_name_update)
    boundary_id = IntProperty(name="Unique ID of This Boundary", default = -1)
    marker = IntProperty(name="Marker Value", default = 1, update=boundary_marker_update)
    status = StringProperty(name="Status")


    def check_boundary_name(self, bnd_name_list):
        """
        Checks for duplicate or illegal boundary name

        @param      self           The object
        @param      bnd_name_list  The bnd name list

        @return     { description_of_the_return_value }
        """

        status = ""
        # Check for duplicate boundary name
        if bnd_name_list.count(self.name) > 1:
            status = "Duplicate boundary: %s" % (self.name)

        # Check for illegal names (Starts with a letter. No special characters)
        bnd_filter = r"(^[A-Za-z]+[0-9A-Za-z_.]*$)"
        m = re.match(bnd_filter, self.name)
        if m is None:
            status = "Boundary name error: %s" % (self.name)

        self.status = status
        return


    def boundary_marker_update(self, context):
        context.active_object['boundaries'][str(self.boundary_id)] = self.marker

    def assign_boundary_faces(self, context):
        # TODO: Remove context from this function...
        obj = context.active_object
        mesh = obj.data

        if mesh.total_face_sel > 0:
            face_set = set()
            with ObjectMode():
                for f in mesh.polygons:
                    if f.select:
                        face_set.add(f.index)

            mats = bpy.data.materials
            bnd_id = self.boundary_id

            bnd_mat = [ mat for mat in mats if mat.gamer.boundary_id == bnd_id ][0]
            act_mat_idx = obj.active_material_index
            bnd_mat_idx = obj.material_slots.find(bnd_mat.name)
            obj.active_material_index = bnd_mat_idx
            bpy.ops.object.material_slot_assign()

            obj.active_material_index = act_mat_idx

            self.set_boundary_faces(context, face_set)
        return {'FINISHED'}


    def repaint_boundary_faces(self, context):
        obj = context.active_object
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')
        self.init_boundary(obj, self.name, self.boundary_id, self.marker)
        self.select_boundary_faces(context)

        bpy.ops.object.mode_set(mode='EDIT')
        mats = bpy.data.materials
        bnd_id = self.boundary_id

        bnd_mat = [ mat for mat in mats if mat.gamer.boundary_id == bnd_id ][0]
        act_mat_idx = obj.active_material_index
        bnd_mat_idx = obj.material_slots.find(bnd_mat.name)
        obj.active_material_index = bnd_mat_idx
        bpy.ops.object.material_slot_assign()
        obj.active_material_index = act_mat_idx


    def remove_boundary_faces(self, context):
        obj = context.active_object
        mesh = obj.data

        if (mesh.total_face_sel > 0):
            face_set = self.get_boundary_faces(context)

            with ObjectMode():
                ml = getMarkerLayer(obj)
                for f in mesh.polygons:
                    if f.select:
                        if f.index in face_set:
                            ml[f.index].value = UNSETID
                        else:   # Remove selection from set to not mess up materials
                            f.select = False
            bnd_unset_mat_idx = obj.material_slots.find('bnd_unset_mat')
            obj.active_material_index = bnd_unset_mat_idx
            bpy.ops.object.material_slot_assign()       # Assign selected material to selected
        return {'FINISHED'}


    def select_boundary_faces(self, context):
        mesh = context.active_object.data
        face_set = self.get_boundary_faces(context)
        with ObjectMode():
            for f in face_set:
                mesh.polygons[f].select = True
        return {'FINISHED'}


    def deselect_boundary_faces(self, context):
        mesh = context.active_object.data
        face_set = self.get_boundary_faces(context)
        with ObjectMode():
            for f in face_set:
                mesh.polygons[f].select = False
        return {'FINISHED'}


    def delete_boundary(self, context):
        """
        Remove boundary data from obj

        @param      self     The object
        @param      context  The context

        @return     { description_of_the_return_value }
        """
        obj = context.active_object
        mesh = obj.data

        with ObjectMode():
            ml = getMarkerLayer(obj)

            # Collect faces and reset marker value
            faces = []
            for i, marker in enumerate(ml):
                if marker.value == self.boundary_id:
                    marker.value = UNSETID
                    faces.append(i)

            idx = 0
            for i in range(0,len(mesh.polygons)):
                if idx <= len(faces)-1:
                    if i == faces[idx]:
                        mesh.polygons[i].select = True
                        idx += 1
                    else:
                        mesh.polygons[i].select = False
                else:
                    mesh.polygons[i].select = False

        mats = bpy.data.materials
        # Assign unset material to members
        bnd_unset_mat_idx = obj.material_slots.find('bnd_unset_mat')
        obj.active_material_index = bnd_unset_mat_idx
        bpy.ops.object.material_slot_assign()       # Assign selected material to selected

        # Remove the boundary object
        obj['boundaries'].pop(str(self.boundary_id))


    def init_boundary(self, obj, bnd_name, bnd_id, bnd_marker):
        """
        Initialize new boundary object

        @param      self
        @param      context     Current Blender context
        @param      bnd_name    String name of the boundary
        @param      bnd_id      String Boundary ID token ('bnd_id_%d')
        @param      bnd_marker  Int boundary marker value
        """
        with ObjectMode():
            ml = getMarkerLayer(obj)

        if not 'boundaries' in obj:
            obj['boundaries'] = dict()      # Initialize boundary object
            obj['boundaries'][str(UNSETID)] = UNSETMARKER      # UNSETID always maps to 0
        obj['boundaries'][str(bnd_id)] = bnd_marker

        mats = bpy.data.materials               # Get list of materials
        if 'bnd_unset_mat' not in mats:
            bnd_unset_mat = bpy.data.materials.new('bnd_unset_mat')
            bnd_unset_mat.use_fake_user = True
            bnd_unset_mat.gamer.boundary_id = UNSETID

        if 'bnd_unset_mat' not in obj.material_slots:
            bpy.ops.object.material_slot_add()  # Add new material slot
            obj.material_slots[0].material = mats['bnd_unset_mat']

        bnd_mat_name = "%s_mat"%(bnd_id)
        if bnd_mat_name not in mats:
            bnd_mat = bpy.data.materials.new(bnd_mat_name)
            # Ensure material is saved even if nothing is allocated to it
            bnd_mat.use_fake_user = True
            bnd_mat.gamer.boundary_id = bnd_id
        else:
            bnd_mat = mats[bnd_mat_name]

        # Add new material to object
        bpy.ops.object.material_slot_add()
        obj.material_slots[-1].material = bnd_mat

        self.boundary_id = bnd_id
        self.marker = bnd_marker
        self.name = bnd_name


    def get_boundary_faces(self, context):
        """
        Given return the set of boundary face indices for this boundary

        @param      self     The object
        @param      context  The context

        @return     The boundary faces.
        """
        face_set = set()
        with ObjectMode():
            ml = getMarkerLayer(context.active_object)
            for i, marker in enumerate(ml):
                if marker.value == self.boundary_id:
                    face_set.add(i)
        return face_set


    def set_boundary_faces(self, context, face_set):
        """
        Set the faces of a given boundary on object, given a set of faces

        @param      self      The object
        @param      context   The context
        @param      face_set  The face set
        """
        obj = context.active_object

        with ObjectMode():
            ml = getMarkerLayer(obj)
            for face in face_set:
                ml[face].value = self.boundary_id


def copyBoundaries(self, fromObject, toObject):
    """
    Copy boundary metadata from boundaries into object
    """
    for bdry in fromObject.gamer.boundary_list:
        toObject.gamer.copyBoundary(toObject, bdry)


class GAMerBoundaryMarkersList(bpy.types.PropertyGroup):
    boundary_list = CollectionProperty(type=GAMerBoundaryMarker, name="Boundary List")
    active_bnd_index = IntProperty(name="Active Boundary Index", default=0)
    get_boundary_info = BoolProperty(
        name="Toggle to enable/disable boundary_info", default=False)


    def get_active_boundary(self):
        """
        @brief      Gets the active boundary object
        @param      self
        @return     GAMerBoundaryMarkerPropertyGroup object of active boundary
        """
        bnd = None
        if len(self.boundary_list) > 0:
            bnd = self.boundary_list[self.active_bnd_index]
        return bnd


    def copyBoundary(self, bdry):
        newBdry = self.boundary_list.add()
        newBdry.init_boundary(bdry.name, bdry.boundary_id, bdry.marker)
        idx = self.boundary_list.find(bnd_name)     # Find index of new boundary
        self.active_bnd_index = idx                 # Set active index to new bdry


    def add_boundary(self, context):
        """
        Add a new boundary to the list of boundaries and set as the active
        boundary
        @param      self     The object
        @param      context  The context
        """
        bnd_id = context.scene.gamer.allocate_boundary_id()
        bnd_name = 'Boundary_%d'%(bnd_id)

        new_bnd = self.boundary_list.add() # Create new GamerBoundaryMarker

        # Name, Boundary ID, Marker
        new_bnd.init_boundary(context.object, bnd_name, bnd_id, bnd_id)

        idx = self.boundary_list.find(bnd_name)     # Find index of new boundary
        self.active_bnd_index = idx                 # Set active index to new bdry


    def repaint_boundaries(self, context):
        for bnd in self.boundary_list:
          bnd.repaint_boundary_faces(context)

        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.object.mode_set(mode='OBJECT')


    def remove_all_boundaries(self, context):
        current_mode = context.object.mode
        for i in range(len(self.boundary_list)):
            # First remove boundary data from mesh:
            bnd = self.boundary_list[0]
            bnd.delete_boundary(context)

            # Now remove the boundary from the object
            self.boundary_list.remove(0)

        self.active_bnd_index = 0 # Restore original mode
        bpy.ops.object.mode_set(mode=current_mode)


    def remove_boundary(self, context):
        current_mode = context.object.mode

        # First remove ID prop boundary data from object:
        bnd = self.get_active_boundary()
        if bnd:
            bnd.delete_boundary(context)

            # Now remove the RNA boundary from the object
            self.boundary_list.remove(self.active_bnd_index)
            self.active_bnd_index -= 1
            if (self.active_bnd_index < 0):
                self.active_bnd_index = 0

        # Restore original mode
        bpy.ops.object.mode_set(mode=current_mode)


    def boundary_name_update(self, context):
        """Performs checks and sorts boundary list after update of boundary names"""
        if self.boundary_list:
            bnd = self.get_active_boundary()
            bnd.check_boundary_name(self.boundary_list.keys())
            self.sort_boundary_list()


    def sort_boundary_list(self):
        """
        Sorts boundary list alphabetically by name
        """

        act_bnd_name = self.get_active_boundary().name

        # Sort the boundary list
        self.inplace_quicksort(self.boundary_list, 0, len(self.boundary_list)-1)

        act_i = self.boundary_list.find(act_bnd_name)
        self.active_bnd_index = act_i


    def inplace_quicksort(self, v, beg, end):  # collection array, int, int
        """
          Sorts a collection array, v, in place.
          Sorts according values in v[i].name
        """
        if ((end - beg) > 0):  # only perform quicksort if we are dealing with > 1 values
            pivot = v[beg].name  # we set the first item as our initial pivot
            i, j = beg, end

            while (j > i):
                while ((v[i].name <= pivot) and (j > i)):
                    i += 1
                while ((v[j].name > pivot) and (j >= i)):
                    j -= 1
                if (j > i):
                    v.move(i, j)
                    v.move(j-1, i)

            if (not beg == j):
                v.move(beg, j)
                v.move(j-1, beg)
            self.inplace_quicksort(v, beg, j-1)
            self.inplace_quicksort(v, j+1, end)


    def draw_panel(self, context, panel):
        layout = panel.layout
        self.draw_layout ( context, layout )

    def draw_layout(self, context, layout):
        active_obj = context.active_object

        if active_obj and (active_obj.type == 'MESH'):
            row = layout.row()
            row.label(text="Defined Boundaries:", icon='FACESEL_HLT')
            row = layout.row()
            col = row.column()
            col.template_list("GAMER_UL_check_boundary", "boundary_list_1",
                          active_obj.gamer, "boundary_list",
                          active_obj.gamer, "active_bnd_index",
                          rows=2)
            col = row.column(align=True)
            col.operator("gamer.add_boundary", icon='ZOOMIN', text="")
            col.operator("gamer.remove_boundary", icon='ZOOMOUT', text="")
            col.operator("gamer.remove_all_boundaries", icon='X', text="")

            # Could have boundary item draw itself in new row here:
            active_bnd = self.get_active_boundary()
            if active_bnd:

                row = layout.row()
                row.prop(active_bnd, "name")

                row = layout.row()
                col = row.column()
                col.label(text="Marker:")
                col = row.column()
                col.prop(active_bnd, "marker", text="")

                mats = bpy.data.materials
                bnd_mat = [ mat for mat in mats if mat.gamer.boundary_id == active_bnd.boundary_id ][0]

            if active_obj.mode == 'EDIT' and active_bnd:
                row = layout.row()
                sub = row.row(align=True)
                sub.operator("gamer.assign_boundary_faces", text="Assign")
                sub.operator("gamer.remove_boundary_faces", text="Remove")

                sub = row.row(align=True)
                sub.operator("gamer.select_boundary_faces", text="Select")
                sub.operator("gamer.deselect_boundary_faces", text="Deselect")
                sub.operator("gamer.select_all_boundary_faces", text="Select All")
