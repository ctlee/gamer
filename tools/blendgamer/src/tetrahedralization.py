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
#
import bpy
from bpy.props import BoolProperty, CollectionProperty, EnumProperty, \
    FloatProperty, FloatVectorProperty, IntProperty, IntVectorProperty, \
    PointerProperty, StringProperty, BoolVectorProperty

import blendgamer.pygamer as g
from blendgamer.util import *

# python imports
import os
import numpy as np


class GAMER_OT_tet_domain_add(bpy.types.Operator):
    bl_idname = "gamer.tet_domain_add"
    bl_label = "Add a Tet Domain"
    bl_description = "Add a new tetrahedralization domain"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.scene.gamer.tet_group.add_tet_domain(context)
        self.report({'INFO'}, "Added a new Tet domain")
        return {'FINISHED'}


class GAMER_OT_tet_domain_remove(bpy.types.Operator):
    bl_idname = "gamer.tet_domain_remove"
    bl_label = "Remove a Tet Domain"
    bl_description = "Remove a tetrahedralization domain"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.scene.gamer.tet_group.remove_active_tet_domain(context)
        self.report({'INFO'}, "Deleted Active Tet Group")
        return {'FINISHED'}


class GAMER_OT_tet_domain_remove_all(bpy.types.Operator):
    bl_idname = "gamer.tet_domain_remove_all"
    bl_label = "Remove a Tet Domain"
    bl_description = "Remove a tetrahedralization domain"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        context.scene.gamer.tet_group.remove_all_tet_domains(context)
        self.report({'INFO'}, "Deleted Active Tet Groups")
        return {'FINISHED'}


class GAMER_OT_tetrahedralize(bpy.types.Operator):
    bl_idname = "gamer.tetrahedralize"
    bl_label = "Tetrahedralize"
    bl_description = ("Tetrahedralize")
    bl_options = {'REGISTER'}

    def execute(self, context):
        context.scene.gamer.tet_group.tetrahedralize(self.report)
        return {'FINISHED'}

    def invoke(self, context, event):
        return self.execute(context)


class GAMerTetDomainPropertyGroup(bpy.types.PropertyGroup):
    # name = StringProperty()  # This is a reminder that "name" is already defined for all subclasses of PropertyGroup
    domain_id = IntProperty(
            name="ID", default=-1,
            description="Domain ID")
    object_name = StringProperty(
            name="ObjName", default="",
            description="Object Name")
    marker = IntProperty(
            name="Marker", default=-1,
            description="Domain Marker Integer")
    is_hole = BoolProperty(
            name="Hole", default=False,
            description="Use this domain as a hole")
    constrain_vol = BoolProperty(
            name="Constrain Volume", default=False,
            description="Constrain Volume")
    vol_constraint = FloatProperty(
            name="Vol Constraint", default=10.0,
            description="Volume Constraint")

    def draw_layout(self, layout):
        row = layout.row()
        col = row.column()
        col.prop(self, "is_hole", text="Use Domain as a Hole")
        if not self.is_hole:
            col = row.column()
            col.prop(self, "marker")
            row = layout.row()
            col = row.column()
            col.prop(self, "constrain_vol")
            if self.constrain_vol:
                col = row.column()
                col.prop(self, "vol_constraint")

    def draw_item_in_row ( self, row ):
        col = row.column()
        col.label(text=str(self.object_name))
        col = row.column()
        col.label(text="Domain ID: " + str(self.domain_id))
        col = row.column()
        if self.is_hole:
            col.label(text="Hole")
        else:
            col.label(text="Domain Marker: " + str(self.marker))

class GAMerTetrahedralizationPropertyGroup(bpy.types.PropertyGroup):
    export_path = StringProperty(
            name="Export Directory",
            description="Path to directory where files will be created",
            default="./", maxlen=1024, subtype='DIR_PATH'
            )
    export_filebase = StringProperty(
            name="Filename",
            description="Base name of the files to export",
            default="gamertetmesh", maxlen=1024, subtype='FILE_NAME'
            )
    domain_list = CollectionProperty(
            type=GAMerTetDomainPropertyGroup, name="Domain List")
    active_domain_index = IntProperty(
            name="Active Domain Index", default=0)
    next_id = IntProperty(
            name="Counter for Unique Domain IDs", default=1)  # Start ID's at 1 to confirm initialization

    show_settings = BoolProperty(
            name="Tetrahedralization Settings", default=False,
            description="Show more detailed settings")

    min_dihedral = FloatProperty(
            name="Min Dihedral", default=10.0,
            description="Minimum Dihedral in Degrees")
    max_aspect_ratio = FloatProperty(
            name="Max Aspect Ratio", default=1.3,
            description="Maximum Aspect Ratio")

    ho_mesh = BoolProperty(
            name="Higher order mesh generation", default=False,
            description="Higher order mesh generation")

    dolfin = BoolProperty(
            name="DOLFIN", default=False,
            description="Generate DOLFIN output")
    # diffpack = BoolProperty(
    #         name="Diffpack", default=False,
    #         description="Generate Diffpack output")
    paraview = BoolProperty(
            name="Paraview", default=False,
            description="Generate Paraview output")

    status = StringProperty ( name="status", default="" )


    def add_tet_domain(self, context):
        # print("Adding a Tet Domain")
        """ Add a new tet domain to the list of tet domains for each selected object """
        #mcell = context.scene.mcell

        # From the list of selected objects, only add MESH objects.
        objs = [obj for obj in context.selected_objects if obj.type == 'MESH']
        if len(objs) > 0:
            for obj in objs:
                # Check by name to see if it's already listed
                current_domain_names = [d.object_name for d in self.domain_list]
                print("Current domains = " + str(current_domain_names))
                if not (obj.name in current_domain_names):
                    new_id = self.allocate_available_id()  # Do this first to check for empty list before adding
                    new_dom = self.domain_list.add()
                    new_dom.domain_id = new_id
                    new_dom.marker = new_id
                    new_dom.object_name = obj.name
                    self.active_domain_index = len(self.domain_list) - 1

    def remove_active_tet_domain(self, context):
        # print("Removing active Tet Domain")
        """ Remove the active tet domain from the list of domains """
        self.domain_list.remove(self.active_domain_index)
        self.active_domain_index -= 1
        if self.active_domain_index < 0:
            self.active_domain_index = 0

    def remove_all_tet_domains(self, context):
        # print("Removing All Tet Domains")
        """ Remove all tet domains from the list of domains """
        while len(self.domain_list) > 0:
            self.domain_list.remove(0)
        self.active_domain_index = 0


    def allocate_available_id(self):
        """ Return a unique domain ID for a new domain """
        # print ("Next ID is " + str(self.next_id))
        if len(self.domain_list) <= 0:
            # Reset the ID to 1 when there are no more molecules
            self.next_id = 1
        self.next_id += 1
        return(self.next_id - 1)

    # def draw_panel(self, context, panel):
    #     layout = panel.layout
    #     self.draw_layout(context, layout)

    def tetrahedralize(self, report):
        print ("######################## Begin Tetrahedralize ########################")

        # filename = self.tet_path
        filename = self.export_path + self.export_filebase
        if not (self.dolfin or self.paraview):
            self.status = "Please select an output format in Tetrahedralization Settings"
            print(self.status)
        else:
            self.status = ""
            mesh_formats = []

            if self.dolfin:
                mesh_formats.append("dolfin")
            if self.paraview:
                mesh_formats.append("paraview")

            # Vector of SurfaceMeshes
            # gmeshes = g.VectorSM()
            gmeshes = list()
            for (obj_name,tet_domain) in [ (d.object_name,d) for d in self.domain_list ]:
                print ( "obj_name = " + obj_name + ", tet_domain = " + str(tet_domain) )

            current_domain_names = [ d.object_name for d in self.domain_list ]
            print ( "Current domains = " + str(current_domain_names) )

            for d in self.domain_list:
                obj = bpy.data.objects[d.object_name]
                gmesh = blenderToGamer(report, obj=obj, map_boundaries=True)
                if not gmesh:
                    print("blenderToGamer returned a gmesh of None")
                else:
                    # # Necessary to prevent garbage collection of gmesh when
                    # # passing into GAMer
                    # gmesh.thisown = 0;

                    # # Collect boundary information
                    # for boundary_name, boundary in zip(boundaries.keys(), boundaries.values()):
                    #     boundary_markers.append((boundary["marker"], boundary_name))

                    print("Mesh %s: num verts: %d numfaces: %d" %(obj_name, gmesh.nVertices, gmesh.nFaces))
                    # Set the domain data on the SurfaceMesh these are the per/domain items as_hole, marker, and volume constraints
                    print("Closed: %d; Marker: %d"%(d.is_hole, d.marker))

                    globalInfo = gmesh.getRoot()
                    globalInfo.ishole = d.is_hole
                    globalInfo.marker = d.marker
                    globalInfo.useVolumeConstraint = d.constrain_vol
                    globalInfo.volumeConstraint = d.vol_constraint
                    # Write surface mesh to file for debug
                    g.writeOFF("surfmesh_%s.off"%(obj_name), gmesh)

                    # Add the mesh
                    gmeshes.append(gmesh)

            # Tetrahedralize mesh
            if len(gmeshes) > 0:

                quality_str = "q%.4f/%.4fO8/7AYVC"%(self.max_aspect_ratio,self.min_dihedral)

                # a%.8f volume constraint..
                quality_str += "o2" if self.ho_mesh else ""

                print("========================================")
                print("TetGen quality string: " + quality_str)
                print("========================================")
                # Do the tetrahedralization
                tetmesh = g.makeTetMesh(gmeshes, quality_str)

                # for i in range(0,5):
                #     print("Laplacian smooth iteration.")
                #     g.smoothMesh(tetmesh)

                # Store mesh to files
                tetmesh_formats  = ["dolfin", "paraview"]
                tetmesh_suffices = [  ".xml", ".vtk"]

                for fmt in mesh_formats:
                    try:
                        suffix = tetmesh_suffices[tetmesh_formats.index(fmt)]
                        print ( "Writing to " + fmt + " file: " + filename + suffix )
                        if fmt == 'dolfin':
                            g.writeDolfin(filename+suffix, tetmesh)
                        if fmt == 'paraview':
                            g.writeVTK(filename+suffix, tetmesh)
                    except Exception as ex:
                        print ( "Error: Unable to write to " + fmt + " file: " + filename + suffix )
                        print ( "   " + str(ex) )

        print ( "######################## End Tetrahedralize ########################" )

classes = [GAMER_OT_tet_domain_add,
           GAMER_OT_tet_domain_remove,
           GAMER_OT_tet_domain_remove_all,
           GAMER_OT_tetrahedralize,
           GAMerTetDomainPropertyGroup,
           GAMerTetrahedralizationPropertyGroup]

def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(make_annotations(cls))

def unregister():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(make_annotations(cls))