#****************************************************************************
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

import blendgamer.report as report
from blendgamer.util import UNSETMARKER

# we use per module class registration/unregistration
def register():
    bpy.utils.register_module(__name__)

def unregister():
    bpy.utils.unregister_module(__name__)


class GAMER_PT_versionerror(bpy.types.Panel):
    bl_label = "GAMer Version Mismatch"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_category = "GAMer"

    @classmethod
    def poll(cls, context):
        return (context.scene is not None) and context.scene.gamer.versionerror != 0

    def draw_header(self, context):
        self.layout.label(text="", icon='CANCEL')

    def draw(self, context):
        layout = self.layout
        layout.alert = True

        if context.scene.gamer.versionerror > 0:
            layout.label(text="Warning the current file was generated with a newer version of GAMer!")
            layout.label(text="We strongly recommend you update the plugin before manipulating this file.")
        elif context.scene.gamer.versionerror < 0:
            layout.label(text="Warning the current file was generated with an older version of GAMer which does not support autoupdate!")
            layout.label(text="Please be careful manipulating this file as data may be lost.")
            layout.label(text="Use the following converters at your own risk!")
            col = layout.column()
            col.operator("gamer.update_to_2_0_1_from_v_0_1")

class GAMER_PT_surfacemesh(bpy.types.Panel):
    bl_label = "Surface Mesh Conditioning"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_category = "GAMer"
    bl_options = {'DEFAULT_CLOSED'}

    # Panel will be drawn if this returns True
    @classmethod
    def poll(cls, context):
        return (context.scene is not None)

    def draw_header(self, context):
        self.layout.label(text="", icon='OUTLINER_DATA_MESH')

    def draw(self, context):
        layout = self.layout
        smprops = context.scene.gamer.surfmesh_procs
        active_obj = context.active_object
        if active_obj and (active_obj.type == 'MESH'):
            box = layout.box()
            col = box.column(align=True)
            if not smprops.advanced_options:
                col.prop(smprops, "advanced_options", icon='TRIA_RIGHT', emboss=False)
            else:
                col.prop(smprops, "advanced_options", icon='TRIA_DOWN', emboss=False)
                col.prop(smprops, "autocorrect_normals")
                col.prop(smprops, "verbose")
                col.prop(smprops, "rings")

            col = layout.column()
            col.label(text="Global mesh operations:")
            col.operator("gamer.normal_smooth", icon='SMOOTHCURVE')
            col.operator("gamer.fill_holes", icon='BORDER_LASSO')
            # col.operator("gamer.refine_mesh", icon='OUTLINER_OB_LATTICE')

            if active_obj.mode == 'EDIT':
                col = layout.column(align=True)
                col.label(text="Local operations (applied to selection only): ")
                rowsub = col.row(align=True)
                rowsub.operator("gamer.coarse_dense",icon='OUTLINER_OB_LATTICE')
                rowsub.prop(smprops, "dense_rate")
                rowsub.prop(smprops, "dense_iter")
                rowsub = col.row(align=True)
                rowsub.operator("gamer.coarse_flat",icon='MOD_TRIANGULATE')
                rowsub.prop(smprops, "flat_rate")
                rowsub.prop(smprops, "flat_iter")

                row = layout.row(align = True)
                row.operator("gamer.smooth", icon='OUTLINER_OB_MESH')
                row.prop(smprops, "smooth_iter")
                row.prop(smprops, "preserve_ridges", expand=True)
            else:
                col = layout.column()
                col.label(text="Change to Edit Mode to enable local improvement options", icon='INFO')
        else:
            layout.label(text="Select a mesh object to use GAMer mesh processing features", icon='HAND')

class GAMER_PT_mesh_quality(bpy.types.Panel):
    bl_label = "Mesh Quality Reporting"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_category = "GAMer"
    bl_options = {'DEFAULT_CLOSED'}

    _type_to_icon = {
    bmesh.types.BMVert: 'VERTEXSEL',
    bmesh.types.BMEdge: 'EDGESEL',
    bmesh.types.BMFace: 'FACESEL',
    }

    # Panel will be drawn if this returns True
    @classmethod
    def poll(cls, context):
        return (context.scene is not None)

    ## This method is from 3D Print Addon
    @staticmethod
    def draw_report(layout, context):
        """Display Reports"""
        info = report.info()
        if info:
            obj = context.edit_object

            layout.label(text="Mesh Stats Report:")
            box = layout.box()
            col = box.column(align=False)
            # box.alert = True
            for i, (text, data) in enumerate(info):
                if obj and data and data[1]:
                    bm_type, bm_array = data
                    col.operator("gamer.meshstats_select_report",
                                 text=text,
                                 icon=GAMER_PT_mesh_quality._type_to_icon[bm_type]).index = i
                else:
                    col.label(text=text)

    def draw(self, context):
        layout = self.layout
        # Mesh quality properties object
        qProps = context.scene.gamer.mesh_quality_properties

        box = layout.box()
        col = box.column()
        col.prop(qProps, "export_path")
        col.prop(qProps, "export_filebase")
        col.operator("gamer.write_quality_info")

        col = layout.column(align=True)
        col.operator("gamer.meshstats_check_solid")
        col.operator("gamer.meshstats_check_degenerate")
        col.operator("gamer.meshstats_check_intersect")

        row = col.row(align=True)
        row.operator("gamer.meshstats_check_wagonwheels")
        row.prop(qProps, "n_wagon_edges")

        row = col.row(align=True)
        row.operator("gamer.meshstats_check_sharp")
        row.prop(qProps, "min_angle")

        col = layout.column()
        col.label(text="Mesh analysis:")
        col.operator("gamer.meshstats_check_all", text="Generate Mesh Report")

        GAMER_PT_mesh_quality.draw_report(layout, context)


        col=layout.column()
        col.label(text="Curvature Calculations:")
        if context.scene.gamer.matplotlib_found:
            row = col.row(align=True)
            row.prop(qProps, "curvePercentile")
            row.prop(qProps, "showplots")
            row.prop(qProps, "saveplots")
            row = col.row(align=True)
            row.prop(qProps, "mixpoint")
            row.prop(qProps, "curveIter")
            row = col.row(align=True)
            row.prop(qProps, "minCurve")
            row.prop(qProps, "maxCurve")


            col.operator("gamer.compute_curvatures")
            col.operator("gamer.mean_curvature")
            col.operator("gamer.gaussian_curvature")
            col.operator("gamer.k1_curvature")
            col.operator("gamer.k2_curvature")

            # col=layout.column()
            # col.label(text="Curvature Calculations:")
            # row = col.row(align=True)
            # row.prop(qProps, "minEnergy")
            # row.prop(qProps, "maxEnergy")
            # col.operator("gamer.helfrich_energy")
        else:
            col.label(text="Curvature Calculations require matplotlib.", icon='LAMP')


class GAMER_PT_boundary_marking(bpy.types.Panel):
    bl_label = "Boundary Marking"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_category = "GAMer"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return (context.scene is not None)

    def draw_header(self, context):
        self.layout.label(text="", icon='TPAINT_HLT')

    def draw(self, context):
        layout = self.layout

        active_obj = context.active_object
        if active_obj and (active_obj.type == 'MESH'):
            row = layout.row()
            row.label("Unmarked marker value = %d"%(UNSETMARKER))
            row = layout.row()
            row.prop(bpy.data.materials['bnd_unset_mat'], 'diffuse_color', text = "Unmarked boundary color")

            row = layout.row()
            row.label(text="Defined Boundaries:", icon='FACESEL_HLT')
            row = layout.row()
            col = row.column()
            col.template_list("GAMER_UL_boundary_list",
                    "GAMer Boundary List",
                    active_obj.gamer, "boundary_list",
                    active_obj.gamer, "active_bnd_index",
                    rows=2,
                    type='DEFAULT'
                )
            col = row.column(align=True)
            col.operator("gamer.add_boundary", icon='ZOOMIN', text="")
            col.operator("gamer.remove_boundary", icon='ZOOMOUT', text="")
            col.operator("gamer.remove_all_boundaries", icon='X', text="")

            # Could have boundary item draw itself in new row here:
            active_bnd = active_obj.gamer.get_active_boundary()
            if active_bnd:
                row = layout.row()
                row.label(text="Set active boundary properties:")

                # Row to update active_bnd name
                row = layout.row()
                row.prop(active_bnd, "boundary_name")
                # Row to update marker value
                row = layout.row()
                row.label(text="Marker:")
                row.prop(active_bnd, "marker", text="") # suppress default txt

                row = layout.row()
                if active_obj.mode == 'OBJECT':
                    row.label(text="Change to Edit Mode to enable boundary assignment", icon='INFO')

            if active_obj.mode == 'EDIT' and active_bnd:
                row = layout.row()
                sub = row.row(align=True)
                sub.operator("gamer.assign_boundary_faces", text="Assign")
                sub.operator("gamer.remove_boundary_faces", text="Remove")

                sub = row.row(align=True)
                sub.operator("gamer.select_boundary_faces", text="Select")
                sub.operator("gamer.deselect_boundary_faces", text="Deselect")

                layout.separator()
                row = layout.row()
                row.operator("gamer.select_all_boundary_faces", text="Select All Marked Boundaries")
        else:
            layout.label(text="Select a mesh object to use GAMer boundary marking features", icon='HAND')

# Object Boundary Panel:
class GAMER_UL_boundary_list(bpy.types.UIList):
    def draw_item(self, context, layout, data, item, icon, active_data,
                  active_propname, index):
        """
        Draw the UI list for boundary markers
        """
        if item.status:
            layout.label(item.boundary_name, icon='ERROR')
        else:
            layout.label(item.boundary_name)

        # Show the color swatch in last section only
        split = layout.split(percentage = 0.5)
        col = split.column()
        col = split.column()
        mats = bpy.data.materials
        bnd_mat = [ mat for mat in mats if mat.gamer.boundary_id == item.boundary_id ][0]
        col.prop(bnd_mat, 'diffuse_color', text='')

class GAMER_PT_tetrahedralization(bpy.types.Panel):
    bl_label = "Tetrahedralization"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'TOOLS'
    bl_category = "GAMer"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return (context.scene is not None)

    def draw_header(self, context):
        self.layout.label(text="", icon='MESH_ICOSPHERE')

    def draw(self, context):
        context.scene.gamer.tet_group.draw_layout(context, self.layout)