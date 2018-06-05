
# ***************************************************************************
# This file is part of the GAMer software.
# Copyright (C) 2016-2017
# by Christopher Lee, John Moody, Rommie Amaro, J. Andrew McCammon,
#    and Michael Holst
# Copyright (C) 2010-2014
# by Zeyun Yu, Michael Holst, Johan Hake, and Tom Bartol

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

bl_info = {
    "name": "GAMer",
    "description": "GAMer: Geometry-preserving Adaptive Mesher",
    "author": "Christopher T. Lee, Zeyun Yu, Michael Holst, Johan Hake, and Tom Bartol",
    "version": (0,1,0),
    "blender": (2, 7, 5),
    "api": 55057,
    "location": "View3D > Add > Mesh",
    "warning": "",
    "wiki_url": "http://www.fetk.org/codes/gamer",
    "tracker_url": "https://github.com/ctlee/gamer/issues",
    "category": "Mesh"}


if "bpy" in locals():
    print("Reloading GAMer")
    import imp
    imp.reload(gamer_gui)
    imp.reload(boundary_markers)
    imp.reload(tetrahedralization)
else:
    print("Importing GAMer")
    from . import gamer_gui
    from . import boundary_markers
    from . import tetrahedralization

# General import
import bpy
import sys
import os

def add_handler ( handler_list, handler_function ):
    """ Only add a handler if it's not already in the list """
    if not (handler_function in handler_list):
        handler_list.append ( handler_function )


def remove_handler ( handler_list, handler_function ):
    """ Only remove a handler if it's in the list """
    if handler_function in handler_list:
        handler_list.remove ( handler_function )


def register():
    print("Registering GAMer...")
    bpy.utils.register_module(__name__)

    bpy.types.Scene.gamer = bpy.props.PointerProperty(
        type=gamer_gui.GAMerPropertyGroup)
    bpy.types.Object.gamer = bpy.props.PointerProperty(
        type=boundary_markers.GAMerBoundaryMarkersListPropertyGroup)
    bpy.types.Material.gamer = bpy.props.PointerProperty(
        type=boundary_markers.GAMerBoundaryMaterialPropertyGroup)

    # Add the load_post handlers
    add_handler ( bpy.app.handlers.load_post, gamer_gui.gamer_load_post )
    add_handler ( bpy.app.handlers.load_post, boundary_markers.boundary_markers_load_post )

    print("GAMer registered")


def unregister():
    remove_handler ( bpy.app.handlers.load_post, boundary_markers.boundary_markers_load_post )
    remove_handler ( bpy.app.handlers.load_post, gamer_gui.gamer_load_post )
    bpy.utils.unregister_module(__name__)

    print("GAMer unregistered")


# for testing
if __name__ == '__main__':
    register()
