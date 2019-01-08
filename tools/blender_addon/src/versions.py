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
from ast import literal_eval

class GAMER_OT_prompt_update(bpy.types.Operator):
    bl_idname = "gamer.prompt_update"
    bl_label = "Warn to update GAMer addon"
    bl_options = {'BLOCKING', 'INTERNAL'}

    def execute(self, context):
        self.report({'WARNING'},
            "Blendfile was generated with a newer version of GAMer.")
        return {'FINISHED'}

def getGamerVersion():
    return sys.modules['gamer_addon'].bl_info.get('version', (-1, -1, -1))

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
    if compare_version(fileVer, currVer) == -1:
        scene.gamer.versionerror = False
        if compare_version(fileVer, (2,0,0)) == 0:
            print("Metadata version is out of date.",
                    "Migrating from v(2, 0, 0) to v%s"%(str(currVer)))
            for obj in bpy.data.objects:
                if obj.type == 'MESH':
                    # Migrate name to boundary_name
                    for bdry in obj.gamer.boundary_list:
                        bdry.boundary_name = bdry.name
                        bdry.name = str(bdry.boundary_id)
                        if 'boundaries' in obj.keys():
                            del obj['boundaries']
            scene.gamer.gamer_version = str(currVer)
    elif compare_version(fileVer, currVer) == 1:
        bpy.ops.gamer.prompt_update()
        scene.gamer.versionerror = True


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
