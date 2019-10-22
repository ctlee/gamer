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
from bpy.props import (BoolProperty, CollectionProperty, EnumProperty,
        FloatProperty, FloatVectorProperty, IntProperty, IntVectorProperty,
        PointerProperty, StringProperty, BoolVectorProperty)
from bpy.app.handlers import persistent

from blendgamer.surfacemesh_ops import SurfaceMeshImprovementProperties
from blendgamer.versions import (checkVersion, getGamerVersion)
from blendgamer.meshstats import MeshQualityReportProperties
from blendgamer.tetrahedralization import GAMerTetrahedralizationPropertyGroup
from blendgamer.markers import GAMerBoundaryMarkersList
from blendgamer.curvatures import GAMerCurvaturesList
from blendgamer.util import UNSETID, make_annotations

import blendgamer.pygamer as pygamer

# python imports
import sys


@persistent
def gamer_load_post(dummy):
    """
    Initialize GAMer addon...
    """
    print('Loading BlendGAMer v%s with PyGAMer %s'%(getGamerVersion(), pygamer.__version__()))
    scene = bpy.context.scene
    scene.gamer.check_for_matplotlib()
    if not scene.gamer.initialized:
        # print('Initializing GAMer Properties')
        scene.gamer.init_properties()
    else:
        # GAMer was previously initialized check version
        checkVersion()
        return

    mats = bpy.data.materials
    if 'bnd_unset_mat' not in mats:
        # if bnd_unset_mat is not defined, then create it
        bnd_unset_mat = bpy.data.materials.new('bnd_unset_mat')
        bnd_unset_mat.use_fake_user = True
        bnd_unset_mat.gamer.boundary_id = UNSETID


class GAMerAddonProperties(bpy.types.PropertyGroup):
    initialized         = BoolProperty(name="GAMer Initialized", default=False)
    matplotlib_found    = BoolProperty(name="Is matplotlib available", default=False)
    gamer_version       = StringProperty(name="GAMer Version", default="0")
    boundary_id_counter = IntProperty(name="GAMer Boundary id Counter")
    versionerror        = IntProperty(name="Version mismatch", default=0)

    surfmesh_procs      = PointerProperty(
                            type=SurfaceMeshImprovementProperties,
                            name="GAMer Surface Mesh Improvement"
                          )

    mesh_quality_properties = PointerProperty(
                            type=MeshQualityReportProperties,
                            name="GAMer Mesh Quality Reporting"
                          )

    tet_group           = PointerProperty(
                            type=GAMerTetrahedralizationPropertyGroup,
                            name="GAMer Tetrahedralization"
                          )

    def allocate_boundary_id ( self ):
        self.boundary_id_counter += 1
        return self.boundary_id_counter

    def init_properties ( self ):
        self.gamer_version = str(getGamerVersion())
        self.boundary_id_counter = 0 # Start counting at 0

        if 'bnd_unset_mat' not in bpy.data.materials:
            bnd_unset_mat = bpy.data.materials.new('bnd_unset_mat')
            bnd_unset_mat.use_fake_user = True
            bnd_unset_mat.gamer.boundary_id = UNSETID
            self.initialized = True

    """"
    Initialize if matplotlib is available
    """
    def check_for_matplotlib(self):
        import importlib.util
        mpl_spec = importlib.util.find_spec("matplotlib")
        self.matplotlib_found = mpl_spec is not None

class GAMerObjectProperties(bpy.types.PropertyGroup):
    markers = PointerProperty(
                type=GAMerBoundaryMarkersList,
                name="Boundary Markers"
              )
    curvatures = PointerProperty(
                    type=GAMerCurvaturesList,
                    name="Curvature Lists"
                 )

classes = [GAMerAddonProperties,
           GAMerObjectProperties]

def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(make_annotations(cls))

def unregister():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(make_annotations(cls))
