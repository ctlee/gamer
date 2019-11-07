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
import sys
import re

from bpy.props import PointerProperty
from ast import literal_eval
from blendgamer.util import *

class GAMER_OT_prompt_update(bpy.types.Operator):
    bl_idname = "gamer.prompt_update"
    bl_label = "Warn to update GAMer addon"
    bl_options = {'BLOCKING', 'INTERNAL'}

    def execute(self, context):
        self.report({'WARNING'},
            "Blendfile was generated with a newer version of GAMer.")
        return {'FINISHED'}

class GAMER_OT_prompt_old_version(bpy.types.Operator):
    bl_idname = "gamer.prompt_old_version"
    bl_label = "Warn that GAMer cannot convert file automatically"
    bl_options = {'BLOCKING', 'INTERNAL'}

    def execute(self, context):
        self.report({'WARNING'},
            "Blendfile was generated with a version that does not support automatic conversion.")
        return {'FINISHED'}

class GAMER_OT_update_to_2_0_1_from_v_0_1(bpy.types.Operator):
    bl_idname = "gamer.update_to_2_0_1_from_v_0_1"
    bl_label = "Update from v0.1 to v2.0.1"
    bl_description = "Update GAMer version to 2.0.1 from v0.1"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        for obj in context.scene.objects:
            # bpy.context.scene.objects.active = obj
            if obj.type == 'MESH':
                print("\n"+"="*30+"\n",obj,"\n"+"="*30)
                if not 'boundaries' in obj:
                    continue
                hidden = False
                if obj.hide:
                    hidden = True
                    obj.hide = False
                obj.gamer.remove_all_boundaries(context)
                for key, bdry in obj['boundaries'].items():
                    print("Migrating boundary: %s"%(key))

                    # First move the material...
                    mat = bpy.data.materials[key+'_mat']
                    newBdryID = context.scene.gamer.boundary_id_counter+1
                    mat.name = materialNamer(newBdryID)
                    mat.gamer.boundary_id = newBdryID
                    mat.use_fake_user = True

                    obj.gamer.add_boundary(context)
                    newBdry = obj.gamer.markers.boundary_list[obj.gamer.active_bnd_index]

                    newBdry.boundary_name = 'NewBoundaryFrom_%s'%(key)
                    newBdry.marker = bdry['marker']

                    # # Deselect all
                    bpy.ops.object.mode_set(mode='EDIT')
                    bpy.ops.mesh.select_all(action='DESELECT')
                    bpy.ops.object.mode_set(mode='OBJECT')

                    # Select faces of interest
                    for faces in bdry['faces'].values():
                        for i in faces:
                            obj.data.polygons[i].select = True
                    bpy.ops.object.mode_set(mode='EDIT')
                    newBdry.assign_boundary_faces(context)
                    bpy.ops.object.mode_set(mode='OBJECT')
                del obj['boundaries']
                if 'id_counter' in obj.gamer:
                    del obj.gamer['id_counter']
                if 'include' in obj.gamer:
                    del obj.gamer['include']
                obj.hide = hidden
        context.scene.gamer.gamer_version = "(2,0,1)"
        checkVersion()
        return {'FINISHED'}


def migrate2_0_1__2_0_5():
    for obj in bpy.data.objects:
        if obj.type == 'MESH':
            # First initialize the data-blocks
            markerProp = obj.gamer.markers
            bl = markerProp.boundary_list

            # Move boundary_list
            if 'boundary_list' in obj['gamer']:
                # Migrate boundary marker info
                obj['gamer']['markers']['boundary_list'] = obj['gamer']['boundary_list']
                del obj['gamer']['boundary_list']
            # Remove active_bnd_index
            if 'active_bnd_index' in obj['gamer']:
                del obj['gamer']['active_bnd_index']
            # Remove include
            if 'include' in obj['gamer']:
                del obj['gamer']['include']

def getGamerVersion():
    return sys.modules['blendgamer'].bl_info.get('version', (-1, -1, -1))

def checkVersion():
    """
    Check the version
    """
    scene = bpy.context.scene
    print("Blendfile contains GAMer v%s metadata"%(scene.gamer.gamer_version))

    fileVer = scene.gamer.gamer_version
    isTupleStr = re.compile('\(.*\)')
    if isTupleStr.match(fileVer):
        fileVer = literal_eval(fileVer)
    else:
        fileVer = tuple(fileVer.split('.'))

    currVer = getGamerVersion()
    scene.gamer.versionerror = compare_version(fileVer, currVer)

    while(scene.gamer.versionerror < 0):
        if scene.gamer.versionerror == -1:
            # Update from 2.0.0 to current
            if compare_version(fileVer, (2,0,0)) == 0:
                newver = (2,0,1)
                print("Metadata version is out of date.",
                        "Migrating from v(2,0,0) to v%s"%(str(newver)))
                for obj in bpy.data.objects:
                    if obj.type == 'MESH':
                        # Migrate name to boundary_name
                        for bdry in obj.gamer.markers.boundary_list:
                            bdry.boundary_name = bdry.name
                            bdry.name = str(bdry.boundary_id)
                            if 'boundaries' in obj.keys():
                                del obj['boundaries']
                scene.gamer.gamer_version = str(newver)
            if compare_version(fileVer, (2,0,1)) >= 0 and compare_version(fileVer, (2,0,5)) < 0 :
                newver = (2,0,5)
                print("Migrating from v%s to v%s"%(str(fileVer), str(newver)))
                migrate2_0_1__2_0_5()
                scene.gamer.gamer_version = str(newver)
            else:
                bpy.ops.gamer.prompt_old_version()
                break
        fileVer = literal_eval(scene.gamer.gamer_version)
        scene.gamer.versionerror = compare_version(fileVer, currVer)

    if scene.gamer.versionerror == 1:
        bpy.ops.gamer.prompt_update()


## VERSION UTILITY FUNCTIONS
def cmp(a, b):
    """
    Compare a and b. Returns -1 if b > a, 1 if a > b, or 0 if a == b
    """
    return (a > b) - (a < b)

def compare_version(v1, v2):
    """
    Compare version tuples

    Return 1:  v1 >  v2
    Return 0:  v1 == v2
    Return -1: v1 <  v2
    """
    return cmp(*zip(*map(lambda x,y:(x or 0, y or 0),
            map(int, v1), map(int, v2))))


classes = [GAMER_OT_prompt_update,
           GAMER_OT_prompt_old_version,
           GAMER_OT_update_to_2_0_1_from_v_0_1]

def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(make_annotations(cls))

def unregister():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(make_annotations(cls))

