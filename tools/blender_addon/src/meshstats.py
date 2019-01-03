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

import bpy
from bpy.types import Operator
from bpy.props import (
    IntProperty,
    FloatProperty,
)
import bmesh

import gamer_addon.report as report
from gamer_addon.util import *

# we use per module class registration/unregistration
def register():
    bpy.utils.register_module(__name__)

def unregister():
    bpy.utils.unregister_module(__name__)


## Following ops are from 3D Print Addon
class MESH_OT_MeshStats_Select_Report(Operator):
    """Select the data associated with this report"""
    bl_idname = "mesh.meshstats_select_report"
    bl_label = "Select Report"
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
        elems = getattr(bm, MESH_OT_MeshStats_Select_Report._type_to_attr[bm_type])[:]

        try:
            for i in bm_array:
                elems[i].select_set(True)
        except:
            # possible arrays are out of sync
            self.report({'WARNING'}, "Report is out of date, re-run check")

        # cool, but in fact annoying
        # bpy.ops.view3d.view_selected(use_all_regions=False)

        return {'FINISHED'}


class MESH_OT_MeshStats_Info_Volume(Operator):
    """Report the volume of the active mesh"""
    bl_idname = "mesh.meshstats_info_volume"
    bl_label = "MeshStats Info Volume"

    @staticmethod
    def main_check(obj,info):
        bm = bmesh_copy_from_object(obj, apply_modifiers=True)
        volume = bm.calc_volume()
        bm.free()
        info.append(("Volume: %s" % clean_float("%.8f" % volume), None))

    def execute(self, context):
        scene = context.scene
        obj = context.active_object

        bm = bmesh_copy_from_object(obj, apply_modifiers=True)
        volume = bm.calc_volume()
        bm.free()

        info = []
        info.append(("Volume: %s" % clean_float("%.8f" % volume), None))

        report.update(*info)
        return {'FINISHED'}


class MESH_OT_MeshStats_Info_Area(Operator):
    """Report the surface area of the active mesh"""
    bl_idname = "mesh.meshstats_info_area"
    bl_label = "MeshStats Info Area"

    @staticmethod
    def main_check(obj, info):
        bm = bmesh_copy_from_object(obj, apply_modifiers=True)
        area = bmesh_calc_area(bm)
        bm.free()
        info.append(("Area: %s" % clean_float("%.8f" % area), None))


    def execute(self, context):
        scene = context.scene
        unit = scene.unit_settings
        obj = context.active_object

        bm = bmesh_copy_from_object(obj, apply_modifiers=True)
        area = bmesh_calc_area(bm)
        bm.free()

        info = []
        info.append(("Area: %s" % clean_float("%.8f" % area), None))

        report.update(*info)
        return {'FINISHED'}


# ---------------
# Geometry Checks

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


class MESH_OT_MeshStats_Check_Solid(Operator):
    """Check for geometry is solid (has valid inside/outside) and correct normals"""
    bl_idname = "mesh.meshstats_check_solid"
    bl_label = "MeshStats Check Solid"

    @staticmethod
    def main_check(obj, info):
        import array

        bm = bmesh_copy_from_object(obj, transform=False, triangulate=False)

        edges_non_manifold = array.array('i', (i for i, ele in enumerate(bm.edges)
                if not ele.is_manifold))
        edges_non_contig = array.array('i', (i for i, ele in enumerate(bm.edges)
                if ele.is_manifold and (not ele.is_contiguous)))

        info.append(("Non Manifold Edge: %d" % len(edges_non_manifold),
                    (bmesh.types.BMEdge, edges_non_manifold)))

        info.append(("Bad Contig. Edges: %d" % len(edges_non_contig),
                    (bmesh.types.BMEdge, edges_non_contig)))

        bm.free()

    def execute(self, context):
        return execute_check(self, context)


class MESH_OT_MeshStats_Check_Intersections(Operator):
    """Check geometry for self intersections"""
    bl_idname = "mesh.meshstats_check_intersect"
    bl_label = "MeshStats Check Intersections"

    @staticmethod
    def main_check(obj, info):
        faces_intersect = bmesh_check_self_intersect_object(obj)
        info.append(("Intersect Face: %d" % len(faces_intersect),
                    (bmesh.types.BMFace, faces_intersect)))

    def execute(self, context):
        return execute_check(self, context)


class MESH_OT_MeshStats_Check_Degenerate(Operator):
    """Check for degenerate geometry that may not print properly """ \
    """(zero area faces, zero length edges)"""
    bl_idname = "mesh.meshstats_check_degenerate"
    bl_label = "MeshStats Check Degenerate"

    @staticmethod
    def main_check(obj, info):
        import array

        scene = bpy.context.scene
        threshold = 0.0001
        # TODO: (0) Update this to store locally in GAMer
        # print_3d = scene.print_3d
        # threshold = print_3d.threshold_zero

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




class MESH_OT_MeshStats_Check_All(Operator):
    """Run all checks"""
    bl_idname = "mesh.meshstats_check_all"
    bl_label = "MeshStats Check All"

    check_cls = (
        MESH_OT_MeshStats_Info_Volume,
        MESH_OT_MeshStats_Info_Area,
        MESH_OT_MeshStats_Check_Solid,
        MESH_OT_MeshStats_Check_Intersections,
        MESH_OT_MeshStats_Check_Degenerate,
        )

    def execute(self, context):
        obj = context.active_object

        info = []

        for cls in self.check_cls:
            cls.main_check(obj, info)

        report.update(*info)

        multiple_obj_warning(self, context)

        return {'FINISHED'}

class MESH_OT_MeshStats_Clean_Isolated(Operator):
    """Cleanup isolated vertices and edges"""
    bl_idname = "mesh.meshstats_clean_isolated"
    bl_label = "MeshStats Clean Isolated "
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        obj = context.active_object
        bm = bmesh_from_object(obj)

        info = []
        change = False

        def face_is_isolated(ele):
            for loop in ele.loops:
                loop_next = loop.link_loop_radial_next
                if loop is not loop_next:
                    return False
            return True

        def edge_is_isolated(ele):
            return ele.is_wire

        def vert_is_isolated(ele):
            return not bool(ele.link_edges)

        # --- face
        elems_remove = [ele for ele in bm.faces if face_is_isolated(ele)]
        remove = bm.faces.remove
        for ele in elems_remove:
            remove(ele)
        change |= bool(elems_remove)
        info.append(("Faces Removed: %d" % len(elems_remove),
                    None))
        del elems_remove
        # --- edge
        elems_remove = [ele for ele in bm.edges if edge_is_isolated(ele)]
        remove = bm.edges.remove
        for ele in elems_remove:
            remove(ele)
        change |= bool(elems_remove)
        info.append(("Edge Removed: %d" % len(elems_remove),
                    None))
        del elems_remove
        # --- vert
        elems_remove = [ele for ele in bm.verts if vert_is_isolated(ele)]
        remove = bm.verts.remove
        for ele in elems_remove:
            remove(ele)
        change |= bool(elems_remove)
        info.append(("Verts Removed: %d" % len(elems_remove),
                    None))
        del elems_remove
        # ---

        report.update(*info)

        if change:
            bmesh_to_object(obj, bm)
            return {'FINISHED'}
        else:
            return {'CANCELLED'}


class MESH_OT_MeshStats_Clean_Non_Manifold(Operator):
    """Cleanup problems, like holes, non-manifold vertices, and inverted normals"""
    bl_idname = "mesh.meshstats_clean_non_manifold"
    bl_label = "MeshStats Clean Non-Manifold and Inverted"
    bl_options = {'REGISTER', 'UNDO'}

    threshold = bpy.props.FloatProperty(
        name="threshold",
        description="Minimum distance between elements to merge",
        default=0.0001,
    )
    sides = bpy.props.IntProperty(
        name="sides",
        description="Number of sides in hole required to fill",
        default=4,
    )

    def execute(self, context):
        self.context = context
        mode_orig = context.mode

        self.setup_environment()
        bm_key_orig = self.elem_count(context)

        self.delete_loose()
        self.remove_doubles(self.threshold)
        self.dissolve_degenerate(self.threshold)

        # may take a while
        self.fix_non_manifold(context, self.sides)

        self.make_normals_consistently_outwards()

        bm_key = self.elem_count(context)

        if mode_orig != 'EDIT_MESH':
            bpy.ops.object.mode_set(mode='OBJECT')

        self.report(
                {'INFO'},
                "Modified Verts:%+d, Edges:%+d, Faces:%+d" %
                (bm_key[0] - bm_key_orig[0],
                 bm_key[1] - bm_key_orig[1],
                 bm_key[2] - bm_key_orig[2],
                 ))

        return {'FINISHED'}

    @staticmethod
    def elem_count(context):
        bm = bmesh.from_edit_mesh(context.edit_object.data)
        return len(bm.verts), len(bm.edges), len(bm.faces)

    @staticmethod
    def setup_environment():
        """set the mode as edit, select mode as vertices, and reveal hidden vertices"""
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_mode(type='VERT')
        bpy.ops.mesh.reveal()

    @staticmethod
    def remove_doubles(threshold):
        """remove duplicate vertices"""
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.remove_doubles(threshold=threshold)

    @staticmethod
    def delete_loose():
        """delete loose vertices/edges/faces"""
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.delete_loose()

    @staticmethod
    def dissolve_degenerate(threshold):
        """dissolve zero area faces and zero length edges"""
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.dissolve_degenerate(threshold=threshold)

    @staticmethod
    def make_normals_consistently_outwards():
        """have all normals face outwards"""
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.normals_make_consistent()

    @classmethod
    def fix_non_manifold(cls, context, sides):
        """naive iterate-until-no-more approach for fixing manifolds"""
        total_non_manifold = cls.count_non_manifold_verts(context)

        if not total_non_manifold:
            return

        bm_states = set()
        bm_key = cls.elem_count(context)
        bm_states.add(bm_key)

        while True:
            cls.fill_non_manifold(sides)

            cls.delete_newly_generated_non_manifold_verts()

            bm_key = cls.elem_count(context)
            if bm_key in bm_states:
                break
            else:
                bm_states.add(bm_key)

    @staticmethod
    def select_non_manifold_verts(
            use_wire=False,
            use_boundary=False,
            use_multi_face=False,
            use_non_contiguous=False,
            use_verts=False,
            ):
        """select non-manifold vertices"""
        bpy.ops.mesh.select_non_manifold(
                extend=False,
                use_wire=use_wire,
                use_boundary=use_boundary,
                use_multi_face=use_multi_face,
                use_non_contiguous=use_non_contiguous,
                use_verts=use_verts,
                )

    @classmethod
    def count_non_manifold_verts(cls, context):
        """return a set of coordinates of non-manifold vertices"""
        cls.select_non_manifold_verts(
                use_wire=True,
                use_boundary=True,
                use_verts=True,
                )

        bm = bmesh.from_edit_mesh(context.edit_object.data)
        return sum((1 for v in bm.verts if v.select))

    @classmethod
    def fill_non_manifold(cls, sides):
        """fill holes and then fill in any remnant non-manifolds"""
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.fill_holes(sides=sides)

        # fill selected edge faces, which could be additional holes
        cls.select_non_manifold_verts(use_boundary=True)
        bpy.ops.mesh.fill()

    @classmethod
    def delete_newly_generated_non_manifold_verts(cls):
        """delete any newly generated vertices from the filling repair"""
        cls.select_non_manifold_verts(use_wire=True, use_verts=True)
        bpy.ops.mesh.delete(type='VERT')
