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
from bpy.props import (
        BoolProperty, CollectionProperty, EnumProperty,
        FloatProperty, FloatVectorProperty, IntProperty, IntVectorProperty,
        PointerProperty, StringProperty, BoolVectorProperty)

from blendgamer.util import *


## OPERATORS
# Object Boundary Marker Operators:
class GAMER_OT_add_boundary(bpy.types.Operator):
    bl_idname       = "gamer.add_boundary"
    bl_label        =  "Add New Boundary"
    bl_description  = "Add a new boundary to an object"
    bl_options      = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.object.gamer.markers.add_boundary(context)
        return {'FINISHED'}


class GAMER_OT_remove_boundary(bpy.types.Operator):
    bl_idname       = "gamer.remove_boundary"
    bl_label        = "Remove Boundary"
    bl_description  = "Remove selected boundary from object"
    bl_options      = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.object.gamer.markers.remove_boundary(context)
        return {'FINISHED'}


class GAMER_OT_remove_all_boundaries(bpy.types.Operator):
    bl_idname       = "gamer.remove_all_boundaries"
    bl_label        = "Remove All Boundaries"
    bl_description  = "Remove all boundaries from object"
    bl_options      = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.object.gamer.markers.remove_all_boundaries(context)
        return {'FINISHED'}


class GAMER_OT_assign_boundary_faces(bpy.types.Operator):
    bl_idname       = "gamer.assign_boundary_faces"
    bl_label        = "Assign Selected Faces To Boundary"
    bl_description  = "Assign selected faces to boundary"
    bl_options      = {'REGISTER', 'UNDO'}

    def execute(self, context):
        bnd = context.object.gamer.markers.get_active_boundary()
        if bnd:
            bnd.assign_boundary_faces(context)
            return {'FINISHED'}
        self.report({'WARNING'},
                "Cannot assign faces, no active boundary selected")
        return {'CANCELLED'}


class GAMER_OT_remove_boundary_faces(bpy.types.Operator):
    bl_idname       = "gamer.remove_boundary_faces"
    bl_label        = "Remove Selected Faces From Boundary"
    bl_description  = "Remove selected faces from boundary"
    bl_options      = {'REGISTER', 'UNDO'}

    def execute(self, context):
        bnd = context.object.gamer.markers.get_active_boundary()
        if bnd:
            bnd.remove_boundary_faces(context)
            return{'FINISHED'}
        self.report({'WARNING'},
                "Cannot remove faces, no active boundary selected")
        return {'CANCELLED'}


class GAMER_OT_select_boundary_faces(bpy.types.Operator):
    bl_idname       = "gamer.select_boundary_faces"
    bl_label        = "Select Faces of Selected Boundary"
    bl_description  = "Select faces of selected boundary"
    bl_options      = {'REGISTER', 'UNDO'}

    def execute(self, context):
        bnd = context.object.gamer.markers.get_active_boundary()
        if bnd:
            bnd.select_boundary_faces(context)
            return {'FINISHED'}
        return {'CANCELLED'}


class GAMER_OT_deselect_boundary_faces(bpy.types.Operator):
    bl_idname       = "gamer.deselect_boundary_faces"
    bl_label        = "Deselect Faces of Selected Boundary"
    bl_description  = "Deselect faces of selected boundary"
    bl_options      = {'REGISTER', 'UNDO'}

    def execute(self, context):
        bnd = context.object.gamer.markers.get_active_boundary()
        if bnd:
            bnd.deselect_boundary_faces(context)
            return {'FINISHED'}
        return {'CANCELLED'}

class GAMER_OT_select_all_boundary_faces(bpy.types.Operator):
    bl_idname       = "gamer.select_all_boundary_faces"
    bl_label        = "Select All Faces of Selected Boundary"
    bl_description  = "Select all faces of selected boundary"
    bl_options      = {'REGISTER', 'UNDO'}

    def execute(self, context):
        for bnd in context.object.gamer.markers.boundary_list:
            bnd.select_boundary_faces(context)
        return {'FINISHED'}

## PROPERTYGROUPS

class GAMerBoundaryMaterial(bpy.types.PropertyGroup):
    """
    @brief      Class for GAMer boundary material property group.
    """
    boundary_id = IntProperty(name="Boundary ID associated with material", default = -1)

class GAMerBoundaryMarker(bpy.types.PropertyGroup):
    """
    @brief      Class for GAMer boundary markers property group.
    """
    boundary_id = IntProperty(
            name = "Boundary ID",
            description = "Unique identifier of this boundary",
            min = -1,
            default = -1
        )
    boundary_name = StringProperty(
            name = "Boundary Name",
            description = "Name of the boundary",
            default = "Boundary"
        )
    marker = IntProperty(
            name = "Marker Value",
            description  = "Marker value to associate with this boundary",
            default = 1
        )
    status = BoolProperty(name="Status", default=False)

    def init_boundary(self, context):
        """
        Initialize new boundary object
        """
        bnd_id = context.scene.gamer.allocate_boundary_id()
        bnd_name = "Boundary_%d"%(bnd_id)

        # Set ID prop key to bnd_id.
        # This should be considered immutable...
        self.name = str(bnd_id)

        obj = context.active_object
        with ObjectMode():
            ml = getMarkerLayer(obj)

        # Get list of materials
        mats = bpy.data.materials
        if 'bnd_unset_mat' not in mats:
            # if bnd_unset_mat is not defined, then create it
            bnd_unset_mat = bpy.data.materials.new('bnd_unset_mat')
            bnd_unset_mat.use_fake_user = True
            bnd_unset_mat.gamer.boundary_id = UNSETID

        if 'bnd_unset_mat' not in obj.material_slots:
            # Add bnd_unset to material_slots...
            bpy.ops.object.material_slot_add()  # Add new material slot
            obj.material_slots[-1].material = mats['bnd_unset_mat']

        bnd_mat_name = materialNamer(bnd_id)
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
        self.marker = bnd_id
        self.boundary_name = bnd_name



    def delete_boundary(self, context):
        """
        Remove boundary data from obj
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

            # Iterate through faces and select accordingly
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

            # Clean up the materials
            mats = bpy.data.materials
            bnd_mat_name = materialNamer(self.boundary_id)
            # First remove from slots
            obj.active_material_index = obj.material_slots.find(bnd_mat_name)
            bpy.ops.object.material_slot_remove()

            # Remove the global materials
            if bnd_mat_name in mats:
                bnd_mat = mats[bnd_mat_name]
                mats.remove(bnd_mat)

            # Assign unset material to members
            bnd_unset_mat_idx = obj.material_slots.find('bnd_unset_mat')
            obj.active_material_index = bnd_unset_mat_idx
            bpy.ops.object.material_slot_assign()


    def assign_boundary_faces(self, context):
        obj = context.active_object
        mesh = obj.data

        if mesh.total_face_sel > 0:
            # Collect list of faces
            face_set = set()
            with ObjectMode():
                for f in mesh.polygons:
                    if f.select:
                        face_set.add(f.index)
            mats = bpy.data.materials
            bnd_id = self.boundary_id

            bnd_mat = getBoundaryMaterial(bnd_id)
            # Cache current active material index
            act_mat_idx = obj.active_material_index

            # Assign material to selected
            bnd_mat_idx = obj.material_slots.find(bnd_mat.name)
            obj.active_material_index = bnd_mat_idx
            bpy.ops.object.material_slot_assign()

            # Reset to cached active material index
            obj.active_material_index = act_mat_idx
            # Set the internal GAMer marker layer values
            self.set_boundary_faces(context, face_set)


    def repaint_boundary_faces(self, context):
        obj = context.active_object
        mats = bpy.data.materials

        bnd_id = self.boundary_id

        if 'bnd_unset_mat' not in obj.material_slots:
            # Add bnd_unset to material_slots...
            bpy.ops.object.material_slot_add()  # Add new material slot
            obj.material_slots[-1].material = mats['bnd_unset_mat']

        bnd_mat_name = materialNamer(bnd_id)
        # Check if material slot is associated
        if bnd_mat_name not in obj.material_slots:
            bpy.ops.object.material_slot_add()  # Add new material slot
            obj.material_slots[-1].material = mats[bnd_mat_name]

        # Deselect all first
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')

        # Select marked faces
        self.select_boundary_faces(context)
        bpy.ops.object.mode_set(mode='EDIT')

        # Repaint boundary
        bnd_mat_idx = obj.material_slots.find(bnd_mat_name)
        obj.active_material_index = bnd_mat_idx
        bpy.ops.object.material_slot_assign()


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


    def select_boundary_faces(self, context):
        obj = context.active_object
        face_set = self.get_boundary_faces(context)

        bm = bmesh_from_object(obj)
        bm.faces.ensure_lookup_table()
        for f in face_set:
            bm.faces[f].select_set(True)

        bmesh_to_object(obj, bm)


    def deselect_boundary_faces(self, context):
        obj = context.active_object

        face_set = self.get_boundary_faces(context)
        bm = bmesh_from_object(obj)
        bm.faces.ensure_lookup_table()
        for f in face_set:
            bm.faces[f].select_set(False)

        bmesh_to_object(obj,bm)

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
        Sets the value of the Boundary Marker Layer
        Set the faces of a given boundary on object, given a set of faces

        @param      context   The context
        @param      face_set  The face set
        """
        obj = context.active_object

        with ObjectMode():
            ml = getMarkerLayer(obj)
            for face in face_set:
                ml[face].value = self.boundary_id


    # def copyBoundaries(self, fromObject, toObject):
    #     """
    #     Copy boundary metadata from boundaries into object
    #     """
    #     for bdry in fromobject.gamer.markers.boundary_list:
    #         toobject.gamer.markers.copyBoundary(toObject, bdry)

    # TODO: (10) Enforce MCell boundary naming conventions
    # def check_boundary_name(self, bnd_name_list):
    #     status = ""
    #     # Check for duplicate boundary name
    #     if bnd_name_list.count(self.name) > 1:
    #         status = "Duplicate boundary: %s" % (self.name)

    #     # Check for illegal names (Starts with a letter. No special characters)
    #     bnd_filter = r"(^[A-Za-z]+[0-9A-Za-z_.]*$)"
    #     m = re.match(bnd_filter, self.name)
    #     if m is None:
    #         status = "Boundary name error: %s" % (self.name)

    #     self.status = status
    #     return



class GAMerBoundaryMarkersList(bpy.types.PropertyGroup):
    boundary_list     = CollectionProperty(
                            type=GAMerBoundaryMarker,
                            name="Boundary List"
                        )
    active_bnd_index  = IntProperty(
                            name="Active Boundary Index",
                            default=0
                        )

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


    def add_boundary(self, context):
        """
        Add a new boundary to the list of boundaries and set as the active
        boundary
        """
        # Create new GamerBoundaryMarker and initialize it
        new_bnd = self.boundary_list.add()
        new_bnd.init_boundary(context)

        # Set active index to new boundary
        idx = self.boundary_list.find(new_bnd.name)
        self.active_bnd_index = idx


    def repaint_boundaries(self, context):
        for bnd in self.boundary_list:
            bnd.repaint_boundary_faces(context)

        # Select all
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.object.mode_set(mode='OBJECT')


    def remove_all_boundaries(self, context):
        for i in range(len(self.boundary_list)):
            # First remove boundary data from mesh:
            bnd = self.boundary_list[0]
            bnd.delete_boundary(context)

            # Now remove the boundary from the object
            self.boundary_list.remove(0)

        self.active_bnd_index = 0 # Restore original mode


    def remove_boundary(self, context):
        # First remove ID prop boundary data from object:
        bnd = self.get_active_boundary()
        if bnd:
            bnd.delete_boundary(context)

            # Now remove the RNA boundary from the object
            self.boundary_list.remove(self.active_bnd_index)
            self.active_bnd_index -= 1
            if (self.active_bnd_index < 0):
                self.active_bnd_index = 0


classes = [GAMerBoundaryMarker,
           GAMerBoundaryMaterial,
           GAMerBoundaryMarkersList,
           GAMER_OT_add_boundary,
           GAMER_OT_remove_boundary,
           GAMER_OT_remove_all_boundaries,
           GAMER_OT_assign_boundary_faces,
           GAMER_OT_remove_boundary_faces,
           GAMER_OT_select_boundary_faces,
           GAMER_OT_deselect_boundary_faces,
           GAMER_OT_select_all_boundary_faces]

def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(make_annotations(cls))

def unregister():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(make_annotations(cls))

