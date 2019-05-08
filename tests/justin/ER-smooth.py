# operations to smoothen ER meshes with very tight geometries

# settings
bpy.context.scene.gamer.surfmesh_procs.smooth_iter = 10
bpy.context.scene.gamer.surfmesh_procs.dense_iter = 5
bpy.context.scene.gamer.surfmesh_procs.dense_rate = 1.2

# set to edit mode
bpy.ops.object.mode_set(mode='EDIT')
bpy.ops.mesh.select_all(action='SELECT')

bpy.ops.mesh.normals_make_consistent(inside=False)


# subdivide
sub_div(2)
smooth_tris(15)
normal_smooth(2)
coarsen(5)
smooth_tris(5)
normal_smooth(1)
smooth_tris(2)





####################





def sub_div(niter):
    for i in range(niter):
        bpy.ops.mesh.subdivide()


def smooth_tris(niter):
    for i in range(niter):
        bpy.ops.gamer.smooth()
        bpy.ops.mesh.select_all(action='SELECT')

def normal_smooth(niter):
    for i in range(niter):
        bpy.ops.gamer.normal_smooth()
        bpy.ops.mesh.select_all(action='SELECT')


def coarsen(niter):
    for i in range(niter):
        bpy.ops.gamer.coarse_dense()
        bpy.ops.mesh.select_all(action='SELECT')


def del_non_manifold(niter):
    """very experimental"""
    for i in range(niter):
        bpy.ops.mesh.select_non_manifold()
        bpy.ops.mesh.delete(type='FACE')
        bpy.ops.mesh.select_more()
        bpy.ops.mesh.select_more()
        bpy.ops.mesh.delete(type='VERT')
        bpy.ops.mesh.select_non_manifold()
        bpy.ops.mesh.delete(type='FACE')
        # fill and triangulate
        bpy.ops.mesh.select_non_manifold()
        bpy.ops.mesh.edge_face_add()
        bpy.ops.mesh.quads_convert_to_tris(quad_method='BEAUTY', ngon_method='BEAUTY')
        bpy.ops.mesh.select_all(action='DESELECT')

