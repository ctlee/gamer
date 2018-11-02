#!/usr/bin/env python

import pygamer as g

# Generate the solvent excluded  molecular surface
mesh = g.ReadPDB_molsurf('2jho.pdb')
# mesh = g.ReadPDB_gauss('2jho.pdb', -0.2, 2.5) # filename, blobbiness, isovalue

# Set outward facing normals
g.compute_orientation(mesh)
g.correctNormals(mesh)

# Get the root node and set some metadata
gInfo = mesh.getRoot()
gInfo.closed = True
gInfo.ishole = True # Don't tetrahedralize this domain
gInfo.marker = -1

# Perform some coarsening and smoothing operations
g.coarse(mesh, 0.016, 1, 0)
g.coarse(mesh, 2.5, 0, 10)
g.smoothMesh(mesh, 15, 165, 5, True)


# Set boundary markers of the mesh to 23
for faceID in mesh.faceIDs():
    faceID.data().marker = 23

# for face in mesh.faces():
#     print(face.marker)

# Generate a surrounding box
box = g.meshCube(5)

# Set box boundary markers to 50
for faceID in box.faceIDs():
    faceID.data().marker = 50

# Get and set box metadata
gInfo = box.getRoot()
gInfo.closed = True
gInfo.ishole = False
gInfo.marker = 5

# Get the radius of the protein mesh and scale the box to be big enough.
radius = g.getRadius(mesh)
g.scale(box, radius*2)


# Generate list of meshes
meshes = g.VectorSM(); # Holder for list of meshes
meshes.push_back(mesh)
meshes.push_back(box)

# Call tetgen
tetmesh = g.MakeTetMesh(meshes, "pq1.4zYAQ")

g.writeVTK('outputTetMesh.vtk', tetmesh)
g.writeDolfin('outputTetMesh.xml', tetmesh)