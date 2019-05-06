# testing gamer vertex smooth

# create cube
bpy.ops.transform.resize(value=(1, 8, 1), constraint_axis=(False, True, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
bpy.ops.mesh.primitive_cube_add(radius=1, view_align=False, enter_editmode=False, location=(0, 0, 0), layers=(True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False))

# scale in y
bpy.ops.transform.resize(value=(1, 8, 1), constraint_axis=(False, True, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)

# scale in x
bpy.ops.transform.resize(value=(0.2, 1, 1), constraint_axis=(True, False, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)

# triangulate
bpy.ops.mesh.quads_convert_to_tris(quad_method='FIXED', ngon_method='CLIP')

# subdivide
nsubd = 5
for i in range(nsubd):
	bpy.ops.mesh.subdivide()

# smooth
nsmooth = 20
for i in range(nsmooth):
	#bpy.ops.mesh.vertices_smooth()
	bpy.ops.gamer.smooth()