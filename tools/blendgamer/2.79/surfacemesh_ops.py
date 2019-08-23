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
from bpy.props import (
        BoolProperty, CollectionProperty, EnumProperty,
        FloatProperty, FloatVectorProperty, IntProperty, IntVectorProperty,
        PointerProperty, StringProperty, BoolVectorProperty)

import blendgamer.pygamer as g
from blendgamer.util import *
from blendgamer.markers import *


# python imports
import os, sys
import numpy as np
import collections

# we use per module class registration/unregistration
def register():
    bpy.utils.register_module(__name__)

def unregister():
    bpy.utils.unregister_module(__name__)


class GAMER_OT_coarse_dense(bpy.types.Operator):
    bl_idname = "gamer.coarse_dense"
    bl_label = "Coarse Dense"
    bl_description = "Decimate selected dense areas of the mesh"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        if context.scene.gamer.surfmesh_procs.coarse_dense(context, self.report):
            self.report({'INFO'}, "GAMer: Coarse Dense complete")
            return {'FINISHED'}
        else:
            return {'CANCELLED'}


class GAMER_OT_coarse_flat(bpy.types.Operator):
    bl_idname = "gamer.coarse_flat"
    bl_label = "Coarse Flat"
    bl_description = "Decimate selected flat areas of the mesh"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        if context.scene.gamer.surfmesh_procs.coarse_flat(context, self.report):
            self.report({'INFO'}, "GAMer: Coarse Flat complete")
            return {'FINISHED'}
        else:
            self.report({"ERROR"}, result)
            return {'CANCELLED'}


class GAMER_OT_smooth(bpy.types.Operator):
    bl_idname = "gamer.smooth"
    bl_label = "Smooth"
    bl_description = "Smooth selected vertices of the mesh"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        if context.scene.gamer.surfmesh_procs.smooth(context, self.report):
            self.report({'INFO'}, "GAMer: Smooth Mesh complete")
            return {'FINISHED'}
        else:
            return {'CANCELLED'}


class GAMER_OT_normal_smooth(bpy.types.Operator):
    bl_idname = "gamer.normal_smooth"
    bl_label = "Normal Smooth"
    bl_description = "Smooth facet normals of selected faces of the mesh"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        if context.scene.gamer.surfmesh_procs.normal_smooth(context, self.report):
            self.report({'INFO'}, "GAMer: Normal Smooth complete")
            return {'FINISHED'}
        else:
            return {'CANCELLED'}

class GAMER_OT_fill_holes(bpy.types.Operator):
    bl_idname = "gamer.fill_holes"
    bl_label = "Fill Holes"
    bl_description = "Triangulate holes"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        if context.scene.gamer.surfmesh_procs.fill_holes(context, self.report):
            self.report({'INFO'}, "GAMer: Fill Holes complete")
            return {'FINISHED'}
        else:
            return {'CANCELLED'}

# class GAMER_OT_refine_mesh(bpy.types.Operator):
#     bl_idname = "gamer.refine_mesh"
#     bl_label = "Quadrisect mesh"
#     bl_description = "Refine the mesh by quadisction"
#     bl_options = {'REGISTER', 'UNDO'}

#     def execute(self, context):
#         if context.scene.gamer.surfmesh_procs.refine_mesh(context, self.report):
#             self.report({'INFO'}, "GAMer: Refine Mesh complete")
#             return {'FINISHED'}
#         else:
#             return {'CANCELLED'}

class SurfaceMeshImprovementProperties(bpy.types.PropertyGroup):
    dense_rate = FloatProperty(
        name="CD_Rate", default=1, min=0.001, max=50.0, precision=4,
        description="The rate for coarsening dense areas")
    dense_iter = IntProperty(
        name="CD_Iter", default=1, min=1, max=15,
        description="The number of iterations for coarsening dense areas")
    flat_rate = FloatProperty(
        name="CF_Rate", default=0.016, min=0.00001, max=0.5, precision=4,
        description="The rate for coarsening flat areas")
    flat_iter = IntProperty(
        name="CF_Iter", default=1, min=1, max=15,
        description="The number of iterations for coarsening flat areas")
    smooth_iter = IntProperty(
        name="S_Iter", default=10, min=1, max=50,
        description="The number of iterations for coarsening dense areas")
    preserve_ridges = BoolProperty(
        name="Preserve ridges", default=False,
        description="Don't flip edges which lie on ridges")
    advanced_options = BoolProperty(name="Advanced options", default=False,
        description="Show additional surface mesh improvement options")
    autocorrect_normals = BoolProperty(name="Autocorrect normals", default=True,
        description="Auto fix inconsistent normals")
    verbose = BoolProperty(name="Verbose", default=False,
        description="Print information to console")
    rings = IntProperty(
        name="LST rings", default=2, min=1, max=5,
        description="The number of neighborhood rings to consider for LST calculation"
        )

    def coarse_dense(self, context, report):
        gmesh = blenderToGamer(report, autocorrect_normals=self.autocorrect_normals)
        if gmesh:
            try:
                gmesh.coarse_dense(rate=self.dense_rate, numiter=self.dense_iter, rings=self.rings)
            except Exception as e:
                report({'ERROR'}, str(e))
                return False
            return gamerToBlender(report, gmesh)
        return False


    def coarse_flat(self, context, report):
        gmesh = blenderToGamer(report, autocorrect_normals=self.autocorrect_normals)
        if gmesh:
            try:
                gmesh.coarse_flat(rate=self.flat_rate, numiter=self.flat_iter, rings=self.rings)
            except Exception as e:
                report({'ERROR'}, str(e))
                return False
            return gamerToBlender(report, gmesh)
        return False


    def smooth(self, context, report):
        gmesh = blenderToGamer(report, autocorrect_normals=self.autocorrect_normals)
        if gmesh:
            try:
                gmesh.smooth(max_iter=self.smooth_iter, preserve_ridges=self.preserve_ridges, rings=self.rings)
            except Exception as e:
                report({'ERROR'}, str(e))
                return False
            return gamerToBlender(report, gmesh)
        return False


    def normal_smooth(self, context, report):
        gmesh = blenderToGamer(report, autocorrect_normals=self.autocorrect_normals)
        if gmesh:
            try:
                gmesh.normalSmooth()
            except Exception as e:
                report({'ERROR'}, str(e))
                return False
            return gamerToBlender(report, gmesh)
        return False


    def fill_holes(self, context, report):
        gmesh = blenderToGamer(report, autocorrect_normals=self.autocorrect_normals)
        if gmesh:
            try:
                gmesh.fillHoles()
            except Exception as e:
                report({'ERROR'}, str(e))
                return False
            return gamerToBlender(report, gmesh)
        return False

    # def refine_mesh(self, context, report):
    #     gmesh = blenderToGamer(report, autocorrect_normals=self.autocorrect_normals)
    #     if gmesh:
    #         try:
    #             g.refine_mesh(gmesh)
    #         except Exception as e:
    #             report({'ERROR'}, str(e))
    #             return False
    #         return gamerToBlender(report, gmesh)
    #     return False
