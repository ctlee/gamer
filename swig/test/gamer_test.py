from gamer import *

# Read a surface mesh and Improve it
sm = SurfaceMesh("single-TT.off")
sm.refine()
print("vertices:", sm.num_vertices, "simplices:", sm.num_faces)
sm.smooth()
sm.coarse_flat()
sm.coarse_dense()
sm.face(0).m = 100
sm.write_off("single-TT_improved.off")

sm = read_PDB_molsurf("1CID.pdb")
print("vertices:", sm.num_vertices, "simplices:", sm.num_faces)

while not sm.smooth(20, 140, 1, True):
    pass

# Improve molecular surface
sm.coarse_flat()
sm.coarse_dense(rate=2.0, numiter=10)
sm.smooth(20, 140, 10, True)
sm.normal_smooth()
sm.remove_unconnected_patches(10)

# Get center information
atom = sm.get_center_radius()

# Generate sphere mesh
outer = SurfaceMesh(5)
ratio = 10
for i in range(outer.num_vertices):
    outer.vertex(i).x = outer.vertex(i).x*atom.radius*ratio+atom.x;
    outer.vertex(i).y = outer.vertex(i).y*atom.radius*ratio+atom.y;
    outer.vertex(i).z = outer.vertex(i).z*atom.radius*ratio+atom.z;

sm.write_off("some_molecule.off")
sm.as_hole = True

gem_mesh = GemMesh([outer, sm])
gem_mesh.write_off("molecular_volmesh.off")
gem_mesh.write_mcsf("molecular_volmesh.m")
gem_mesh.write_dolfin("molecular_volmesh.xml")
gem_mesh.write_diffpack("molecular_volmesh.grid", [])
