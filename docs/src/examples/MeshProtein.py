#!/usr/bin/env python

import pygamer as g

# Generate the solvent excluded  molecular surface
mesh = g.readPDB_molsurf('2jho.pdb')
# mesh = g.readPDB_gauss('2jho.pdb', -0.2, 2.5) # filename, blobbiness, isovalue

# Set outward facing normals
mesh.compute_orientation()
mesh.correctNormals()

# Get the root node and set some metadata
gInfo = mesh.getRoot()
gInfo.ishole = True # Don't tetrahedralize this domain
gInfo.marker = -1

# Perform some coarsening and smoothing operations
mesh.coarse(0.016, 1, 0)
mesh.coarse(2.5, 0, 10)

# Set selection of all vertices to True so smooth will operate on them.
for v in mesh.vertexIDs:
    v.data().selected = True
# mesh.smooth(3, True, True)

# Set boundary markers of the mesh to 23
for faceID in mesh.faceIDs:
    faceID.data().marker = 23

# Generate a surrounding box
box = g.surfacemesh.cube(5)

# Set box boundary markers to 50
for faceID in box.faceIDs:
    faceID.data().marker = 50

# Get and set box metadata
gInfo = box.getRoot()
gInfo.ishole = False
gInfo.marker = 5

# Get the radius of the protein mesh and scale the box to be big enough.
center, radius = mesh.getCenterRadius()
box.scale(radius*2)

g.writeOFF('2jho.off', mesh)
g.writeOFF('box.off', box)

# Construct a list of meshes to pass into TetGen
meshes = [mesh, box]

# Call tetgen to tetrahedralize
tetmesh = g.makeTetMesh(meshes, "q1.3/10O8/7AVC")

g.writeVTK('outputTetMesh.vtk', tetmesh)
g.writeDolfin('outputTetMesh.xml', tetmesh)