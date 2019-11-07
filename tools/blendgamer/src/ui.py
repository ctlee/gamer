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
from blendgamer.util import UNSETMARKER, make_annotations
import blendgamer.pygamer as pygamer


if bpy.app.version < (2,80,0):
    REGION = "TOOLS"
    BULB_ICON = 'LAMP'
    ADD_ICON = 'ZOOMIN'
    REMOVE_ICON = 'ZOOMOUT'
    LOVE_ICON = 'COLOR_RED'
else:
    REGION = "UI"
    BULB_ICON = 'LIGHT'
    ADD_ICON = 'ADD'
    REMOVE_ICON = 'REMOVE'
    LOVE_ICON = 'FUND'


class GAMER_PT_versionerror(bpy.types.Panel):
    bl_label = "GAMer Version Mismatch"
    bl_space_type = 'VIEW_3D'
    bl_region_type = REGION
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
    bl_region_type = REGION
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
            col.operator("gamer.normal_smooth")
            col.operator("gamer.fill_holes")

            if active_obj.mode == 'EDIT':
                col = layout.column()
                col.label(text="Local operations (applied to selection only): ")

                col = layout.column(align=True)
                col.operator("gamer.coarse_dense")
                rowsub = col.row(align=True)
                rowsub.prop(smprops, "dense_rate")
                rowsub.prop(smprops, "dense_iter")


                col = layout.column(align=True)
                col.operator("gamer.coarse_flat")
                rowsub = col.row(align=True)
                rowsub.prop(smprops, "flat_rate")
                rowsub.prop(smprops, "flat_iter")


                col = layout.column(align=True)
                col.operator("gamer.smooth")
                rowsub = col.row(align = True)
                rowsub.prop(smprops, "smooth_iter")
                rowsub.prop(smprops, "preserve_ridges", expand=True)
            else:
                col = layout.column()
                col.label(text="Change to Edit Mode to enable local improvement options", icon='INFO')
        else:
            layout.label(text="Select a mesh object to use GAMer mesh processing features", icon='HAND')

# class GAMER_PT_advanced_options(bpy.types.Panel):
#     bl_label = "Advanced Options"
#     bl_parent_id = 'GAMER_PT_surfacemesh'
#     bl_space_type = 'VIEW_3D'
#     bl_region_type = REGION
#     bl_category = "GAMer"
#     bl_options = {'DEFAULT_CLOSED'}


#     def draw(self, context):
#         layout = self.layout

#         smprops = context.scene.gamer.surfmesh_procs
#         active_obj = context.active_object

#         col = layout.column(align=True)
#         col.prop(smprops, "autocorrect_normals")
#         col.prop(smprops, "verbose")
#         col.prop(smprops, "rings")


class GAMER_PT_mesh_quality(bpy.types.Panel):
    bl_label = "Mesh Quality Reporting"
    bl_space_type = 'VIEW_3D'
    bl_region_type = REGION
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
        col.label(text="Curvature Estimation:")

        if context.scene.gamer.matplotlib_found:
            obj = context.object
            if obj.type == 'MESH':
                curveProp = obj.gamer.curvatures
                row = col.row(align=True)
                row.operator("gamer.compute_curvatures")
                row.prop(curveProp, "algorithm", text="")

                row = layout.row()
                row.label(text="Computed Curvatures:", icon='FACESEL')
                row = layout.row()
                col = row.column()
                col.template_list("GAMER_UL_curvature_list",
                        "GAMer Curvature List",
                        curveProp, "curvature_list",
                        curveProp, "active_index",
                        rows=2,
                        type='DEFAULT'
                    )
                col = row.column(align=True)
                col.operator("gamer.remove_curvature", icon=REMOVE_ICON, text="")
                col.operator("gamer.remove_all_curvatures", icon='X', text="")

                crv = curveProp.get_active_index()
                if (crv != None):
                    col = layout.column()
                    row = col.row()
                    row.label(text="Set active curvature properties:")

                    row = col.row()
                    row.prop(crv, "colormap")
                    row.prop(crv, "limitsArePercentiles")

                    col = layout.column(align=True)
                    row = col.row(align=True)
                    row.prop(crv, "mixpoint")
                    row.prop(crv, "curveIter")
                    row = col.row(align=True)
                    row.prop(crv, "minCurve")
                    row.prop(crv, "maxCurve")

                    row = layout.row()
                    row.label(text="Plot Settings:")
                    row = layout.row()
                    row.prop(curveProp, "showplots")
                    row.prop(curveProp, "saveplots")

                    row = layout.row()
                    row.operator("gamer.plot_curvature")
                    row.operator("gamer.plot_all_curvatures")
                    # row = layout.row()
                    # row.operator("gamer.plot_differences")
            else:
                col.label(text="Select a mesh object to enable estimation of curvatures", icon='LIGHT')
        else:
            col.label(text="Curvature estimations require matplotlib.", icon='LIGHT')

# Object Boundary Panel:
class GAMER_UL_curvature_list(bpy.types.UIList):
    def draw_item(self, context, layout, data, item, icon, active_data,
                  active_propname, index):
        """
        Draw the UI list for boundary markers
        """
        row = layout.row()
        row.label(text=layout.enum_item_description(item, 'algorithm', item.algorithm))
        # row.prop_enum(item, 'algorithm', item.algorithm, text=layout.enum_item_description(item, 'algorithm', item.algorithm))
        row.prop_enum(item, 'curvatureType', item.curvatureType)


class GAMER_PT_boundary_marking(bpy.types.Panel):
    bl_label = "Boundary Marking"
    bl_space_type = 'VIEW_3D'
    bl_region_type = REGION
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
            markerProp = active_obj.gamer.markers

            row = layout.row()
            row.label(text="Unmarked marker value = %d"%(UNSETMARKER))
            row = layout.row()
            row.prop(bpy.data.materials['bnd_unset_mat'], 'diffuse_color', text = "Unmarked boundary color")

            row = layout.row()
            row.label(text="Defined Boundaries:", icon='FACESEL')
            row = layout.row()
            col = row.column()
            col.template_list("GAMER_UL_boundary_list",
                    "GAMer Boundary List",
                    markerProp, "boundary_list",
                    markerProp, "active_bnd_index",
                    rows=2,
                    type='DEFAULT'
                )
            col = row.column(align=True)
            col.operator("gamer.add_boundary", icon=ADD_ICON, text="")
            col.operator("gamer.remove_boundary", icon=REMOVE_ICON, text="")
            col.operator("gamer.remove_all_boundaries", icon='X', text="")

            # Could have boundary item draw itself in new row here:
            active_bnd = markerProp.get_active_boundary()
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
            layout.label(text=item.boundary_name, icon='ERROR')
        else:
            layout.label(text=item.boundary_name)

        # Show the color swatch in last section only
        if bpy.app.version < (2,80,0):
            split = layout.split(percentage = 0.5, align = True)
        else:
            split = layout.split(factor = 0.5, align = True)
        col = split.column()
        col = split.column()
        mats = bpy.data.materials
        bnd_mat = [ mat for mat in mats if mat.gamer.boundary_id == item.boundary_id ][0]
        col.prop(bnd_mat, 'diffuse_color', text='')


class GAMER_UL_domain(bpy.types.UIList):
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

        item.draw_item_in_row(layout.row())


class GAMER_PT_tetrahedralization(bpy.types.Panel):
    bl_label = "Tetrahedralization"
    bl_space_type = 'VIEW_3D'
    bl_region_type = REGION
    bl_category = "GAMer"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return (context.scene is not None)

    def draw_header(self, context):
        self.layout.label(text="", icon='MESH_ICOSPHERE')

    def draw(self, context):
        layout = self.layout

        tetprops = context.scene.gamer.tet_group

        row = layout.row()
        row.label(text="Domains")

        row = layout.row()
        col = row.column()

        col.template_list("GAMER_UL_domain", "",
                          tetprops, "domain_list",
                          tetprops, "active_domain_index",
                          rows=2)

        col = row.column(align=True)
        col.operator("gamer.tet_domain_add", icon=ADD_ICON, text="")
        col.operator("gamer.tet_domain_remove", icon=REMOVE_ICON, text="")
        col.operator("gamer.tet_domain_remove_all", icon='X', text="")

        if len(tetprops.domain_list) > 0:
            domain = tetprops.domain_list[tetprops.active_domain_index]

            # row = layout.row()
            # row.label ( "Active Index = " + str ( self.active_domain_index ) + ", ID = " + str ( domain.domain_id ) )

            domain.draw_layout ( layout )

            box = layout.box()
            row = box.row(align=True)
            row.alignment = 'LEFT'
            if not tetprops.show_settings:
                row.prop(tetprops, "show_settings", icon='TRIA_RIGHT', emboss=False)
            else:
                row.prop(tetprops, "show_settings", icon='TRIA_DOWN', emboss=False)

                row = box.row()
                row.prop(tetprops, "export_path")
                row = box.row()
                row.prop(tetprops, "export_filebase")

                row = box.row()
                col = row.column()
                col.prop(tetprops, "min_dihedral")
                col = row.column()
                col.prop(tetprops, "max_aspect_ratio")

                # row = box.row()
                # row.prop ( self, "ho_mesh" )

                row = box.row()
                row.label(text="Output Formats:")

                row = box.row()
                sbox = row.box()

                row = sbox.row()
                col = row.column()
                col.prop(tetprops, "dolfin")
                col = row.column()
                col.prop(tetprops, "paraview")

            row = layout.row()
            icon = 'PROP_OFF'
            if tetprops.dolfin or tetprops.paraview:
                icon = 'PROP_ON'
            row.operator("gamer.tetrahedralize", text="Tetrahedralize", icon=icon)
            if len(tetprops.status) > 0:
                row = layout.row()
                row.label(text=tetprops.status, icon="ERROR")


class GAMER_PT_version(bpy.types.Panel):
    bl_label = "BlendGAMer Info"
    bl_space_type = 'VIEW_3D'
    bl_region_type = REGION
    bl_category = "GAMer"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return (context.scene is not None)

    def draw_header(self, context):
        self.layout.label(text="", icon=LOVE_ICON)

    def draw(self, context):
        layout = self.layout
        layout.alert = True

        layout.label(text="BlendGAMer %s"%(pygamer.__version__()))
        layout.operator("wm.url_open", text="How to Acknowledge").url = "https://gamer.readthedocs.io/en/latest/#acknowleding-the-use-of-gamer-in-your-work"

classes = [GAMER_PT_versionerror,
           GAMER_PT_surfacemesh,
           # GAMER_PT_advanced_options,
           GAMER_PT_mesh_quality,
           GAMER_UL_curvature_list,
           GAMER_PT_boundary_marking,
           GAMER_UL_boundary_list,
           GAMER_UL_domain,
           GAMER_PT_tetrahedralization,
           GAMER_PT_version]


def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(make_annotations(cls))

def unregister():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(make_annotations(cls))

