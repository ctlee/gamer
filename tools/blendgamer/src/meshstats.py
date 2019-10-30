# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

import math
import bpy
from bpy.props import (
        BoolProperty, CollectionProperty, EnumProperty,
        FloatProperty, FloatVectorProperty, IntProperty, IntVectorProperty,
        PointerProperty, StringProperty, BoolVectorProperty)
import bmesh

import blendgamer.pygamer as g

import importlib.util
mpl_spec = importlib.util.find_spec("matplotlib")
mpl_found = mpl_spec is not None

if mpl_found:
    from blendgamer.colormap import dataToVertexColor

from blendgamer.colormap_enums import colormap_enums

import blendgamer.report as report
from blendgamer.util import *


## Following ops are from 3D Print Addon
class GAMER_OT_MeshStats_Select_Report(bpy.types.Operator):
    """Select the data associated with this report"""
    bl_idname  = "gamer.meshstats_select_report"
    bl_label   = "Select Report"
    bl_options = {'INTERNAL'}

    index = IntProperty()

    _type_to_mode = {
        bmesh.types.BMVert: 'VERT',
        bmesh.types.BMEdge: 'EDGE',
        bmesh.types.BMFace: 'FACE',
        }

    _type_to_attr = {
        bmesh.types.BMVert: "verts",
        bmesh.types.BMEdge: "edges",
        bmesh.types.BMFace: "faces",
        }

    def execute(self, context):
        obj = context.edit_object
        info = report.info()
        text, data = info[self.index]
        bm_type, bm_array = data

        bpy.ops.mesh.reveal()
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.mesh.select_mode(type=self._type_to_mode[bm_type])

        bm = bmesh.from_edit_mesh(obj.data)
        elems = getattr(bm, GAMER_OT_MeshStats_Select_Report._type_to_attr[bm_type])[:]

        try:
            for i in bm_array:
                elems[i].select_set(True)
        except:
            # possible arrays are out of sync
            self.report({'WARNING'}, "Report is out of date, re-run check")
        return {'FINISHED'}


# Helper method to get object and run main_check
def execute_check(self, context):
    obj = context.active_object
    info = []
    self.main_check(obj, info)
    report.update(*info)
    multiple_obj_warning(self, context)
    return {'FINISHED'}

def multiple_obj_warning(self, context):
    if len(context.selected_objects) > 1:
        self.report({"INFO"}, "Multiple selected objects. Only the active one will be evaluated")

class GAMER_OT_MeshStats_Info_Volume(bpy.types.Operator):
    """Report the volume of the active mesh"""
    bl_idname  = "gamer.meshstats_info_volume"
    bl_label   = "MeshStats Info Volume"

    @staticmethod
    def main_check(obj,info):
        bm = bmesh_copy_from_object(obj)
        volume = 0.0
        for face in bm.faces:
            if len(face.loops) != 3:
                info.append(("Cannot compute volume for non triangulated object.", None))
                return
            tv0 = face.loops[0].vert.co
            tv1 = face.loops[1].vert.co
            tv2 = face.loops[2].vert.co
            x0 = tv0.x
            y0 = tv0.y
            z0 = tv0.z
            x1 = tv1.x
            y1 = tv1.y
            z1 = tv1.z
            x2 = tv2.x
            y2 = tv2.y
            z2 = tv2.z
            det = x0*(y1*z2-y2*z1)+x1*(y2*z0-y0*z2)+x2*(y0*z1-y1*z0)
            volume = volume + det
        bm.free()
        volume = volume/6.0
        info.append(("Volume: %s" % clean_float("%.8f" % volume), None))

    def execute(self, context):
        return execute_check(self, context)


class GAMER_OT_MeshStats_Info_Area(bpy.types.Operator):
    """Report the surface area of the active mesh"""
    bl_idname   = "gamer.meshstats_info_area"
    bl_label    = "MeshStats Info Area"

    @staticmethod
    def main_check(obj, info):
        bm = bmesh_copy_from_object(obj)
        area = bmesh_calc_area(bm)
        bm.free()
        info.append(("Area: %s" % clean_float("%.8f" % area), None))

    def execute(self, context):
        return execute_check(self, context)

class GAMER_OT_MeshStats_Check_Solid(bpy.types.Operator):
    """Check for geometry is solid (has valid inside/outside) and correct normals"""
    bl_idname   = "gamer.meshstats_check_solid"
    bl_label    = "MeshStats Check Solid"
    bl_description = "Check for non-manifolds and inconsistent normals"

    @staticmethod
    def main_check(obj, info):
        import array

        bm = bmesh_copy_from_object(obj, transform=False, triangulate=False)

        edges_non_manifold = array.array('i', (i for i, ele in enumerate(bm.edges)
                if not ele.is_manifold))
        edges_non_contig = array.array('i', (i for i, ele in enumerate(bm.edges)
                if ele.is_manifold and (not ele.is_contiguous)))

        verts_non_manifold = array.array('i', (i for i, ele in enumerate(bm.verts)
                if not ele.is_manifold))

        info.append(("Non Manifold Edge: %d" % len(edges_non_manifold),
                    (bmesh.types.BMEdge, edges_non_manifold)))

        info.append(("Bad Contig. Edges: %d" % len(edges_non_contig),
                    (bmesh.types.BMEdge, edges_non_contig)))

        info.append(("Non Manifold Vertices: %d" % len(verts_non_manifold),
            (bmesh.types.BMVert, verts_non_manifold)))
        bm.free()

    def execute(self, context):
        return execute_check(self, context)


class GAMER_OT_MeshStats_Check_Intersections(bpy.types.Operator):
    """Check geometry for self intersections"""
    bl_idname   = "gamer.meshstats_check_intersect"
    bl_label    = "MeshStats Check Intersections"

    @staticmethod
    def main_check(obj, info):
        faces_intersect = bmesh_check_self_intersect_object(obj)
        info.append(("Intersect Face: %d" % len(faces_intersect),
                    (bmesh.types.BMFace, faces_intersect)))

    def execute(self, context):
        return execute_check(self, context)


class GAMER_OT_MeshStats_Check_Degenerate(bpy.types.Operator):
    """Check for degenerate geometry that may not print properly """ \
    """(zero area faces, zero length edges)"""
    bl_idname    = "gamer.meshstats_check_degenerate"
    bl_label     = "Check Degenerate Faces and Edges"
    bl_description  = "Check for zero length/area edges and faces"

    @staticmethod
    def main_check(obj, info):
        import array

        threshold = 0
        bm = bmesh_copy_from_object(obj, transform=False, triangulate=False)

        faces_zero = array.array('i', (i for i, ele in enumerate(bm.faces) if ele.calc_area() <= threshold))
        edges_zero = array.array('i', (i for i, ele in enumerate(bm.edges) if ele.calc_length() <= threshold))

        info.append(("Zero Area Faces: %d" % len(faces_zero),
                    (bmesh.types.BMFace, faces_zero)))

        info.append(("Zero Len. Edges: %d" % len(edges_zero),
                    (bmesh.types.BMEdge, edges_zero)))

        bm.free()

    def execute(self, context):
        return execute_check(self, context)


class GAMER_OT_MeshStats_Check_Wagonwheels(bpy.types.Operator):
    bl_idname    = "gamer.meshstats_check_wagonwheels"
    bl_label     = "Check for wagon wheels"
    bl_description = "Check for vertices connected to many edges"

    @staticmethod
    def main_check(obj, info):
        import array

        n_wagon_edges = bpy.context.scene.gamer.mesh_quality_properties.n_wagon_edges
        bm = bmesh_copy_from_object(obj, transform=False, triangulate=False)

        wagon_edges = array.array(
                'i',
                (i for i, ele in enumerate(bm.verts)
                    if len(ele.link_edges) >= n_wagon_edges)
            )
        info.append((
                "Number of Wagonwheels: %d"%len(wagon_edges),
                (bmesh.types.BMVert, wagon_edges)
            ))
        bm.free()

    def execute(self, context):
        return execute_check(self,context)


class GAMER_OT_MeshStats_Check_Sharp(bpy.types.Operator):
    bl_idname   = "gamer.meshstats_check_sharp"
    bl_label    = "Check for small angles"
    bl_description = "Check for faces with small angles"

    @staticmethod
    def main_check(obj, info):
        import array

        min_angle = bpy.context.scene.gamer.mesh_quality_properties.min_angle
        bm = bmesh_copy_from_object(obj, transform=False, triangulate=False)

        # TODO: (10) make this into a fancy list comprehension...
        sharp_list = []
        for i, face in enumerate(bm.faces):
            for loop in face.loops:
                if loop.calc_angle()*180/math.pi <= min_angle:
                    sharp_list.append(i)

        sharp = array.array('i', sharp_list)

        info.append((
                "Sharp faces: %d"%len(sharp),
                (bmesh.types.BMFace, sharp)
            ))
        bm.free()

    def execute(self, context):
        return execute_check(self,context)


class GAMER_OT_MeshStats_Betti_Numbers(bpy.types.Operator):
    bl_idname   = "gamer.meshstats_betti"
    bl_label    = "Report Betti numbers"
    bl_description = "Compute the first three Betti numbers"

    @staticmethod
    def main_check(obj, info):
        gmesh = blenderToGamer(report, autocorrect_normals=False)
        if gmesh:
            # k, orientable, manifold = gmesh.check_orientation()
            # holes = gmesh.nEdges - gmesh.nVertices + k

            valid, k, h, v = gmesh.getBettiNumbers()

            info.append(("Euler Characteristic: %d"%(gmesh.nVertices-gmesh.nEdges+gmesh.nFaces), None))

            if valid:
                info.append(("B0 Connected Components: %d"%(k), None))
                info.append(("B1 Holes: %d"%(h), None))
                info.append(("B2 Voids: %d"%(v), None))
            else:
                info.append(("B0 Connected Components: %d"%(k), None))
                info.append(("Higher order Betti numbers undetermined", None))

    def execute(self, context):
        return execute_check(self,context)


class GAMER_OT_MeshStats_Check_All(bpy.types.Operator):
    """Run all checks"""
    bl_idname = "gamer.meshstats_check_all"
    bl_label  = "MeshStats Check All"
    bl_description = "Check all"

    def execute(self, context):
        obj = context.active_object
        info = []

        check_classes = (
            GAMER_OT_MeshStats_Info_Volume,
            GAMER_OT_MeshStats_Info_Area,
            GAMER_OT_MeshStats_Check_Wagonwheels,
            GAMER_OT_MeshStats_Check_Sharp,
            GAMER_OT_MeshStats_Check_Solid,
            GAMER_OT_MeshStats_Check_Intersections,
            GAMER_OT_MeshStats_Check_Degenerate,
            GAMER_OT_MeshStats_Betti_Numbers,
        )

        for cls in check_classes:
            cls.main_check(obj, info)

        report.update(*info)
        multiple_obj_warning(self, context)
        return {'FINISHED'}


class GAMER_OT_write_quality_info(bpy.types.Operator):
    bl_idname      = "gamer.write_quality_info"
    bl_label       = "Print mesh quality info to files"
    bl_description = "Dump quality info to files"
    bl_options     = {'REGISTER'}

    def execute(self, context):
        mqp = bpy.context.scene.gamer.mesh_quality_properties
        for obj in context.selected_objects:
            if obj.type == 'MESH':
                fname = mqp.export_path + mqp.export_filebase + "_" + obj.name
                print("Dumping quality info of mesh %s to file %s"%(obj.name, fname))
                gmesh = blenderToGamer(self.report, obj=obj)
                if gmesh:
                    g.printQualityInfo(fname, gmesh)
        return {'FINISHED'}


class MeshQualityReportProperties(bpy.types.PropertyGroup):
    n_wagon_edges = IntProperty(
        name="N Edges", default=8, min=1,
        description="The number of incident edges to a vertex to be selected")
    export_path = StringProperty(
            name="Export Directory",
            description="Path to directory where files will be created",
            default="./", maxlen=1024, subtype='DIR_PATH'
            )
    export_filebase = StringProperty(
            name="Filename",
            description="Base name of the files to export",
            default="meshquality", maxlen=1024, subtype='FILE_NAME'
            )
    min_angle = IntProperty(
        name="Angle Threshold", default=15, min=0, max=180,
        description="Select faces with angles less than this criteria")

classes = [
    GAMER_OT_MeshStats_Select_Report,
    GAMER_OT_MeshStats_Info_Volume,
    GAMER_OT_MeshStats_Info_Area,
    GAMER_OT_MeshStats_Check_Solid,
    GAMER_OT_MeshStats_Check_Intersections,
    GAMER_OT_MeshStats_Check_Degenerate,
    GAMER_OT_MeshStats_Check_Wagonwheels,
    GAMER_OT_MeshStats_Check_Sharp,
    GAMER_OT_MeshStats_Betti_Numbers,
    GAMER_OT_MeshStats_Check_All,
    GAMER_OT_write_quality_info,
    MeshQualityReportProperties]

def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(make_annotations(cls))

def unregister():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(make_annotations(cls))