# ***************************************************************************
# This file is part of the GAMer software.
# Copyright (C) 2016-2021
# by Christopher T. Lee and contributors
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, see <http:#www.gnu.org/licenses/>
# or write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# ***************************************************************************

import bpy
from bpy.props import (
    BoolProperty,
    CollectionProperty,
    EnumProperty,
    FloatProperty,
    FloatVectorProperty,
    IntProperty,
    IntVectorProperty,
    PointerProperty,
    StringProperty,
    BoolVectorProperty,
)
from bpy.app.handlers import persistent

from blendgamer.surfacemesh_ops import SurfaceMeshImprovementProperties
from blendgamer.versions import checkVersion, getGamerVersion
from blendgamer.meshstats import MeshQualityReportProperties
from blendgamer.tetrahedralization import GAMerTetrahedralizationPropertyGroup
from blendgamer.markers import GAMerBoundaryMarkersList
from blendgamer.curvatures import GAMerCurvaturesList
from blendgamer.util import UNSETID, make_annotations, get_material_by_bnd_id, getBndUnsetMat

import blendgamer.pygamer as pygamer

# python imports
import sys


@persistent
def gamer_load_post(dummy):
    """
    Initialize the GAMer addon.

    Checks that the addon properties are correctly initialized. Also that the
    plugin version matches the metadata layout stored in a blendfile. This
    function also sets a flag noting whether matplotlib is installed or not.
    """
    print(
        "Loading BlendGAMer v%s with PyGAMer %s"
        % (getGamerVersion(), pygamer.__version__())
    )
    scene = bpy.context.scene
    scene.gamer.check_for_matplotlib()
    if not scene.gamer.initialized:
        # print('Initializing GAMer Properties')
        scene.gamer.init_properties()
    else:
        # GAMer was previously initialized check version
        checkVersion()
        return

    bnd_unset_mat = getBndUnsetMat()


class GAMerAddonProperties(bpy.types.PropertyGroup):
    """
    Property group to store GAMer addon metadata
    """

    initialized = BoolProperty(name="GAMer Initialized", default=False)
    matplotlib_found = BoolProperty(name="Is matplotlib available", default=False)
    gamer_version = StringProperty(name="GAMer Version", default="0")
    boundary_id_counter = IntProperty(name="GAMer Boundary id Counter")
    versionerror = IntProperty(name="Version mismatch", default=0)

    surfmesh_improvement_properties = PointerProperty(
        type=SurfaceMeshImprovementProperties,
        name="GAMer Surface Mesh Improvement Properties",
    )

    mesh_quality_properties = PointerProperty(
        type=MeshQualityReportProperties, name="GAMer Mesh Quality Reporting"
    )

    tet_group = PointerProperty(
        type=GAMerTetrahedralizationPropertyGroup, name="GAMer Tetrahedralization"
    )

    def allocate_boundary_id(self):
        """Allocate the next available boundary ID.

        Returns
        -------
        int
            The boundary ID.
        """
        self.boundary_id_counter += 1
        return self.boundary_id_counter

    def init_properties(self):
        """Initialize BlendGAMer addon properties"""
        self.gamer_version = str(getGamerVersion())
        self.boundary_id_counter = 0  # Start counting at 0

        if "bnd_unset_mat" not in bpy.data.materials:
            bnd_unset_mat = bpy.data.materials.new("bnd_unset_mat")
            bnd_unset_mat.use_fake_user = True
            bnd_unset_mat.gamer.boundary_id = UNSETID
            self.initialized = True

    def check_for_matplotlib(self):
        """Check if matplotlib is available and set an internal flag."""
        import importlib.util

        mpl_spec = importlib.util.find_spec("matplotlib")
        self.matplotlib_found = mpl_spec is not None


class GAMerObjectProperties(bpy.types.PropertyGroup):
    """
    PropertyGroup of properties to link into Blender Objects
    """

    markers = PointerProperty(type=GAMerBoundaryMarkersList, name="Boundary Markers")
    curvatures = PointerProperty(type=GAMerCurvaturesList, name="Curvature Lists")


classes = [GAMerAddonProperties, GAMerObjectProperties]


def register():
    for cls in classes:
        bpy.utils.register_class(make_annotations(cls))


def unregister():
    for cls in reversed(classes):
        bpy.utils.unregister_class(make_annotations(cls))
