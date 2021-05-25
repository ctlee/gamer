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
# ***************************************************************************

import bpy
import bmesh
from bpy.props import (
    BoolProperty,
    CollectionProperty,
    EnumProperty,
    FloatProperty,
    FloatVectorProperty,
    IntProperty,
    IntVectorProperty,
    PointerProperty,
    StringProperty,
    BoolVectorProperty,
)

from blendgamer.util import *


# OPERATORS
# Object Boundary Marker Operators:
class GAMER_OT_add_boundary(bpy.types.Operator):
    bl_idname = "gamer.add_boundary"
    bl_label = "Add New Boundary"
    bl_description = "Add a new boundary to an object"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        context.object.gamer.markers.add_boundary(context)
        return {"FINISHED"}


class GAMER_OT_remove_boundary(bpy.types.Operator):
    bl_idname = "gamer.remove_boundary"
    bl_label = "Remove Boundary"
    bl_description = "Remove selected boundary from object"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        context.object.gamer.markers.remove_boundary(context)
        return {"FINISHED"}


class GAMER_OT_remove_all_boundaries(bpy.types.Operator):
    bl_idname = "gamer.remove_all_boundaries"
    bl_label = "Remove All Boundaries"
    bl_description = "Remove all boundaries from object"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        context.object.gamer.markers.remove_all_boundaries(context)
        return {"FINISHED"}


class GAMER_OT_assign_boundary_faces(bpy.types.Operator):
    bl_idname = "gamer.assign_boundary_faces"
    bl_label = "Assign Selected Faces To Boundary"
    bl_description = "Assign selected faces to boundary"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        bnd = context.object.gamer.markers.get_active_boundary()
        if bnd:
            bnd.assign_boundary_faces(context)
            return {"FINISHED"}
        self.report({"WARNING"}, "Cannot assign faces, no active boundary selected")
        return {"CANCELLED"}


class GAMER_OT_remove_boundary_faces(bpy.types.Operator):
    bl_idname = "gamer.remove_boundary_faces"
    bl_label = "Remove Selected Faces From Boundary"
    bl_description = "Remove selected faces from boundary"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        bnd = context.object.gamer.markers.get_active_boundary()
        if bnd:
            bnd.remove_boundary_faces(context)
            return {"FINISHED"}
        self.report({"WARNING"}, "Cannot remove faces, no active boundary selected")
        return {"CANCELLED"}


class GAMER_OT_select_boundary_faces(bpy.types.Operator):
    bl_idname = "gamer.select_boundary_faces"
    bl_label = "Select Faces of Selected Boundary"
    bl_description = "Select faces of selected boundary"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        bnd = context.object.gamer.markers.get_active_boundary()
        if bnd:
            bnd.select_boundary_faces(context)
            return {"FINISHED"}
        return {"CANCELLED"}


class GAMER_OT_deselect_boundary_faces(bpy.types.Operator):
    bl_idname = "gamer.deselect_boundary_faces"
    bl_label = "Deselect Faces of Selected Boundary"
    bl_description = "Deselect faces of selected boundary"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        bnd = context.object.gamer.markers.get_active_boundary()
        if bnd:
            bnd.deselect_boundary_faces(context)
            return {"FINISHED"}
        return {"CANCELLED"}


class GAMER_OT_select_all_boundary_faces(bpy.types.Operator):
    bl_idname = "gamer.select_all_boundary_faces"
    bl_label = "Select All Faces of Selected Boundary"
    bl_description = "Select all faces of selected boundary"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        for bnd in context.object.gamer.markers.boundary_list:
            bnd.select_boundary_faces(context)
        return {"FINISHED"}


class GAMER_OT_deselect_all_boundary_faces(bpy.types.Operator):
    bl_idname = "gamer.deselect_all_boundary_faces"
    bl_label = "Deselect all marked faces"
    bl_description = "Deselect all faces of selected boundary"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        for bnd in context.object.gamer.markers.boundary_list:
            bnd.deselect_boundary_faces(context)
        return {"FINISHED"}


class GAMerBoundaryMaterial(bpy.types.PropertyGroup):
    """
    Class for GAMer boundary material property group.
    """

    boundary_id = IntProperty(name="Boundary ID associated with material", default=-1)


class GAMerBoundaryMarker(bpy.types.PropertyGroup):
    """
    Class for GAMer boundary markers property group.
    """

    boundary_id = IntProperty(
        name="Boundary ID",
        description="Unique identifier of this boundary",
        min=-1,
        default=-1,
    )
    boundary_name = StringProperty(
        name="Boundary Name", description="Name of the boundary", default="Boundary"
    )
    marker = IntProperty(
        name="Marker Value",
        description="Marker value to associate with this boundary",
        default=1,
    )
    status = BoolProperty(name="Status", default=False)

    def init_boundary(self, context):
        """
        Initialize new boundary object
        """
        bnd_id = context.scene.gamer.allocate_boundary_id()
        bnd_name = "Boundary_%d" % (bnd_id)

        # Set ID prop key to bnd_id.
        # This should be considered immutable...
        self.name = str(bnd_id)

        obj = context.active_object

        with BMeshContext(obj) as bm:
            ml = getBMeshMarkerLayer

        # Get list of materials
        mats = bpy.data.materials
        bnd_unset_mat = getBndUnsetMat()
        # Ensure bnd_unset_mat is in object.material_slots
        if bnd_unset_mat.name not in obj.material_slots:
            # Add bnd_unset to material_slots...
            obj.data.materials.append(bnd_unset_mat)

        # Create a new boundary material. It shouldn't exist already
        # unless GAMerAddonProperties has been tampered with...
        bnd_mat_name = materialNamer(bnd_id)
        bnd_mat = bpy.data.materials.new(bnd_mat_name)
        # Ensure material is saved even if nothing is allocated to it
        bnd_mat.use_fake_user = True
        bnd_mat.gamer.boundary_id = bnd_id

        # Add new material to object material slots
        obj.data.materials.append(bnd_mat)

        self.boundary_id = bnd_id
        self.marker = bnd_id
        self.boundary_name = bnd_name

    def delete_boundary(self, context):
        """
        Remove boundary data from obj
        """
        obj = context.active_object
        mesh = obj.data

        # Clean up the materials
        mats = bpy.data.materials
        bnd_mat = get_material_by_bnd_id(self.boundary_id)

        with ObjectMode():
            # Material slots can only be removed in object mode!
            objmats = obj.data.materials
            # First remove all instances from slots
            while bnd_mat.name in objmats:
                idx = objmats.find(bnd_mat.name)
                if bpy.app.version < (2, 81, 0):
                    objmats.pop(index=idx, update_data=True)
                else:
                    objmats.pop(index=idx)

        # Remove the global material
        mats.remove(bnd_mat)

        # Assign unset material to members
        bnd_unset_mat = getBndUnsetMat()
        bnd_unset_mat_idx = obj.material_slots.find(bnd_unset_mat.name)
        # Link the material to the object if it's somehow missing
        if bnd_unset_mat_idx == -1:
            obj.data.materials.append(bnd_unset_mat)
            bnd_unset_mat_idx = obj.material_slots.find(bnd_unset_mat.name)

        with BMeshContext(obj) as bm:
            ml = getBMeshMarkerLayer(bm)

            for face in bm.faces:
                if face[ml] == self.boundary_id:
                    face[ml] = UNSETID
                    face.material_index = bnd_unset_mat_idx

    def assign_boundary_faces(self, context):
        """Assign boundary marker to selected faces

        Args:
            context (TYPE): Blender context
        """
        obj = context.active_object
        mesh = obj.data

        # Material to associate
        bnd_mat = get_material_by_bnd_id(self.boundary_id)
        matID = obj.material_slots.find(bnd_mat.name)
        # Link the material to the object if it's somehow missing
        if matID == -1:
            obj.data.materials.append(bnd_mat)
            matID = obj.material_slots.find(bnd_mat.name)

        if mesh.total_face_sel > 0:
            with BMeshContext(obj) as bm:
                ml = getBMeshMarkerLayer(bm)
                for face in bm.faces:
                    # Apply boundary marker label if selected
                    if face.select == True:
                        face[ml] = self.boundary_id
                        face.material_index = matID

    def repaint_boundary_faces(self, context):
        obj = context.active_object
        mats = bpy.data.materials

        # Add unset boundary to materials
        bnd_mat = getBndUnsetMat()
        matID = obj.material_slots.find(bnd_mat.name)
        # Link the material to the object if it's somehow missing
        if matID == -1:
            obj.data.materials.append(bnd_mat)
            matID = obj.material_slots.find(bnd_mat.name)

        # Material to associate
        bnd_mat = get_material_by_bnd_id(self.boundary_id)
        matID = obj.material_slots.find(bnd_mat.name)
        # Link the material to the object if it's somehow missing
        if matID == -1:
            obj.data.materials.append(bnd_mat)
            matID = obj.material_slots.find(bnd_mat.name)

        with BMeshContext(obj) as bm:
            ml = getBMeshMarkerLayer(bm)
            for face in bm.faces:
                # Apply boundary marker material
                if face[ml] == self.boundary_id:
                    face.material_index = matID

    def remove_boundary_faces(self, context):
        obj = context.active_object
        mesh = obj.data

        # Material to associate
        bnd_mat = getBndUnsetMat()
        matID = obj.material_slots.find(bnd_mat.name)
        # Link the material to the object if it's somehow missing
        if matID == -1:
            obj.data.materials.append(bnd_mat)
            matID = obj.material_slots.find(bnd_mat.name)

        if mesh.total_face_sel > 0:
            with BMeshContext(obj) as bm:
                ml = getBMeshMarkerLayer(bm)
                for face in bm.faces:
                    if face.select and face[ml] == self.boundary_id:
                        face[ml] = UNSETID
                        face.material_index = matID

    def select_boundary_faces(self, context):
        obj = context.active_object

        with BMeshContext(obj) as bm:
            ml = getBMeshMarkerLayer(bm)

            for face in bm.faces:
                if face[ml] == self.boundary_id:
                    face.select_set(True)

    def deselect_boundary_faces(self, context):
        obj = context.active_object

        with BMeshContext(obj) as bm:
            ml = getBMeshMarkerLayer(bm)

            for face in bm.faces:
                if face[ml] == self.boundary_id:
                    face.select_set(False)
            # bm.select_flush_mode()

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
    boundary_list = CollectionProperty(type=GAMerBoundaryMarker, name="Boundary List")
    active_bnd_index = IntProperty(name="Active Boundary Index", default=0)

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

    def remove_all_boundaries(self, context):
        for i in range(len(self.boundary_list)):
            # First remove boundary data from mesh:
            bnd = self.boundary_list[0]
            bnd.delete_boundary(context)

            # Now remove the boundary from the object
            self.boundary_list.remove(0)

        self.active_bnd_index = 0  # Restore original mode

    def remove_boundary(self, context):
        # First remove ID prop boundary data from object:
        bnd = self.get_active_boundary()
        if bnd:
            bnd.delete_boundary(context)

            # Now remove the RNA boundary from the object
            self.boundary_list.remove(self.active_bnd_index)
            self.active_bnd_index -= 1
            if self.active_bnd_index < 0:
                self.active_bnd_index = 0


classes = [
    GAMerBoundaryMarker,
    GAMerBoundaryMaterial,
    GAMerBoundaryMarkersList,
    GAMER_OT_add_boundary,
    GAMER_OT_remove_boundary,
    GAMER_OT_remove_all_boundaries,
    GAMER_OT_assign_boundary_faces,
    GAMER_OT_remove_boundary_faces,
    GAMER_OT_select_boundary_faces,
    GAMER_OT_deselect_boundary_faces,
    GAMER_OT_select_all_boundary_faces,
    GAMER_OT_deselect_all_boundary_faces,
]


def register():
    from bpy.utils import register_class

    for cls in classes:
        register_class(make_annotations(cls))


def unregister():
    from bpy.utils import unregister_class

    for cls in reversed(classes):
        unregister_class(make_annotations(cls))
