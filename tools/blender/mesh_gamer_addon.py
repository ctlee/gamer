# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 08:37:50 2012

@author: Ludovic Autin <autin@scripps.edu>
"""
__url__ = ["gamer",]

bl_info = {
    "name": "GAMer",
    "description": """""",
    "author": "",
    "version": (0,0,0),
    "blender": (2, 6, 1),
    "api": 31236,
    "location": "View3D > Add > Mesh",
    "warning": '', # used for warning icon and text in addons panel
    "wiki_url": ""\
        "Scripts/My_Script",
    "tracker_url": ""\
        "func=detail&aid=<number>",
    "category": "Mesh"}

# -------------------------------------------------------------------------- 
# ***** BEGIN GPL LICENSE BLOCK ***** 
# 
# This program is free software; you can redistribute it and/or 
# modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation; either version 2 
# of the License, or (at your option) any later version. 
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License 
# along with this program; if not, write to the Free Software Foundation, 
# Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. 
# 
# ***** END GPL LICENCE BLOCK ***** 
# -------------------------------------------------------------------------- 

# General import
import bpy
import sys
import os

bpath = bpy.app.binary_path

# GamerUI should be in the python path too
# also need to add path to upy

#PLUGIN_ID = 1027431 this for c4d
import upy
upy.setUIClass()

plugTypeClass, opType = upy.getPluginClass(plug="command")#= operator in blender

from gamer.upy_gui import GamerUI

# print (bpy.types.Operator)
# print ((plugTypeClass))

class gamer_plugin(plugTypeClass):
    def setgui(self,dname):
        self.gui = GamerUI(title = "GAMer")
        self.gui.setup()
        self.hasGui = True
        self.gui.display()
        
    def resetgui(self,dname):
        self.gui = GamerUI(title = "GAMer")
        self.gui.setup()
        self.gui.display()

# problem maya use the class
gameraddon = gamer_plugin(name="GAMer", pluginId=222222,
                          tooltip="GAMer mesh improvments, mesh annotations "\
                          "and volumetric mesh generation",
                          hasGui=True)

# gamer.setIcon(image_name="gamer.tif")
if "__res__" in locals() :
    gameraddon.register(gamer_plugin, Object=gameraddon,
                        menuadd={"head":None, "mtmesh":None}, res=__res__)
else :
    gameraddon.register(gamer_plugin, Object=gameraddon,
                        menuadd={"head":None, "mtmesh":None})

def register():
    print (__name__)

def unregister():
    pass

#if __name__ == "__main__":
#    epmv_plugin = epmv_Dialog()
#    register()

