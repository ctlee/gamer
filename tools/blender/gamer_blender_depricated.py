#!BPY
"""
Name: 'GAMer mesh improvments'
Blender: 249b
Group: 'Mesh'
Tooltip: 'GAMer mesh improvments'
"""

######################################################
# Importing modules
######################################################
 
import numpy as np
import Blender
import BPyMessages
from Blender import NMesh
from Blender.BGL import *
from Blender.Draw import *

import math, sys, gamer, re

# GAMer parameters
GAMER_SUFFIX = "gi"
MAX_MESH_LENGTH_NAME = 16
NUM_TICK_FIGURES = 2
gparams = dict()
gparams["coarse_dense"] = dict(rate=Create(1.6), numiter=Create(1))
gparams["coarse_flat"] = dict(rate=Create(0.016))
gparams["smooth"] = dict(max_min_angle = Create(15),
                         min_max_angle = Create(150),
                         max_iter = Create(6),
                         preserve_ridges=Create(0))
gparams["normal_smooth"] = dict()
overwrite_mesh = 0
material_button = Create(0)
marker_value = Create(1)
TETMESH_OPTIONS = ["FEtk"]
TETMESH_SUFFIX = [".m"]
tetmesh_value = Create(0)

# Events
EVENT_NOEVENT = 1
EVENT_SMOOTH = 2
EVENT_COARSE_FLAT = 3
EVENT_COARSE_DENSE = 4
EVENT_NORMAL_SMOOTH = 5
EVENT_MATERIAL_SET = 6
EVENT_UPDATE_MIN_MAX_ANGLE = 7
EVENT_UPDATE_MARKER_VALUE = 8
EVENT_TETRAHEDRALIZE = 9
EVENT_TOGGLE_PRESERVE_RIDGES = 10
EVENT_EIGENVALUES = 11
EVENT_REPAIR_MESH = 12
EVENT_EXIT = 20
EVENT_TOGGLE_OVERWRITE = 30
GAMER_TEXT_X = 0.
GAMER_TEXT_Y = 0.
GAMER_MSG = ""

MARKER_DEFAULT = 1
YSTART = 300
BUTTON_HEIGHT = 20
BUTTON_WIDTH = 100
XTEXT = 10
XBUTTON = 10
NUMBER_HEIGHT = 20
NUMBER_WIDTH = 110
DY = max(NUMBER_HEIGHT, BUTTON_HEIGHT) + 0
DX = 0
NUMBER_WIDTH_AND_DX = NUMBER_WIDTH+DX
BUTTON_WIDTH_AND_DX = BUTTON_WIDTH+DX

markers = dict()

def save_to_registry():
    "Save values to registry"
    d = {}
    d["tetmesh_value"] = tetmesh_value.val
    Blender.Registry.SetKey('gamer_params', d, True)
    # Save material markers
    d = {}
    for mtrl in markers:
        d[mtrl.getName()] = markers[mtrl]
    Blender.Registry.SetKey('material_markers', d, True)
        

def load_from_registry():
    "Get values from registry"
    global tetmesh_value, markers
    d = Blender.Registry.GetKey('gamer_params', True)
    if d:
        try:
            tetmesh_value.val = d["tetmesh_value"]
        except: save_to_registry() # If data is not valid, rewrite it.

    d = Blender.Registry.GetKey('material_markers', True)
    if d:
        mtrl_names, mtrls = get_materials()
        try:
            for name in d:
                if name in mtrl_names:
                    mtrl = mtrls[mtrl_names.index(name)]
                    markers[mtrl] = d[name]
        except: save_to_registry() # If data is not valid, rewrite it.

def get_materials():
    "Get all defined materials"
    scn = Blender.Scene.GetCurrent()
    obj = scn.objects.active
    if not obj or obj.type != 'Mesh':
        return [], []
    
    bmesh = obj.getData(mesh=1)
    mtrls = []
    mtrl_names = []
    for i, mat in enumerate(bmesh.materials):
        if mat is not None:
            #print i, mat.getName()
            mtrls.append(mat)
            mtrl_names.append(mat.getName())
    return mtrl_names, mtrls

def select_material():
    global markers, marker_value
    mtrl_names, mtrls = get_materials()
    if material_button.val < len(mtrls):
        marker_value.val = markers[mtrls[material_button.val]]
    Redraw(1)
    
def update_marker_value():
    global markers 
    mtrl_names, mtrls = get_materials()
    if markers:
        markers[mtrls[material_button.val]] = marker_value.val
    
# Load any saved markers
load_from_registry()

# Update the marker value from the loaded markers
mtrl_names, mtrls = get_materials()
if markers and mtrls:
    marker_value.val = markers[mtrls[material_button.val]]
update_marker_value()

######################################################
# GUI drawing
######################################################
def draw():
    global gparams, GAMER_TEXT_X, GAMER_TEXT_Y
 
    ypos = YSTART

        ########## Titles
    glClear(GL_COLOR_BUFFER_BIT)
    glColor3f(0, 0, 0)
    glRasterPos2d(XTEXT, ypos)
    Text("GAMer Improvements", "large")
    ypos -= int(DY*1.5)
    # FIXME: Inplace imnprovement does not work...
    #draw_overwrite(XBUTTON, ypos)
    #ypos -= DY
    draw_repair_mesh(XBUTTON, ypos)
    ypos -= DY
    draw_coarse_dense(XBUTTON, ypos)
    ypos -= DY
    draw_coarse_flat(XBUTTON, ypos)
    ypos -= DY
    draw_smooth(XBUTTON, ypos)
    ypos -= DY
    Button("Normal Smooth", EVENT_NORMAL_SMOOTH, XBUTTON, ypos,
           BUTTON_WIDTH, BUTTON_HEIGHT, "Smooth the normals of the mesh")
    ypos -= DY
    Button("Eigenvalues", EVENT_EIGENVALUES, XBUTTON, ypos,
           BUTTON_WIDTH, BUTTON_HEIGHT, "Print out the eigenvalues of the local structure tensor")
    ypos -= DY
    glRasterPos2d(XTEXT, ypos)
    Text("Assign boundary markers to materials", "large")
    ypos -= int(DY*1.5)
    draw_materials(XBUTTON, ypos)
    ypos -= DY
    draw_tetrahedralize(XBUTTON, ypos)
    #Text(GAMER_MSG)
    ypos -= DY
    ypos -= DY
    Button("Exit", EVENT_EXIT, XBUTTON, ypos, BUTTON_WIDTH, BUTTON_HEIGHT)
    ypos -= DY

def exit_event():
    save_to_registry()
    Exit()

def event(evt, val):
    if evt == QKEY and not val:
        exit_event()
        
def bevent(evt):
    ######### Manages GUI events
    if evt in event_dict:
        event_dict[evt]()

def draw_overwrite(xpos, ypos):
    Toggle("overwrite", EVENT_TOGGLE_OVERWRITE, xpos, ypos, BUTTON_WIDTH,
           BUTTON_HEIGHT, overwrite_mesh, "Overwrite present mesh")
    
def draw_repair_mesh(xpos, ypos):
    Button("Repair mesh", EVENT_REPAIR_MESH, xpos, ypos, BUTTON_WIDTH,
           BUTTON_HEIGHT, "Prepair mesh for GAMer")

def draw_coarse_dense(xpos, ypos):
    global gparams
    Button("Coarse Dense", EVENT_COARSE_DENSE, xpos, ypos, BUTTON_WIDTH,
           BUTTON_HEIGHT, "Coarse dense areas of the mesh")
    gparams["coarse_dense"]["rate"] =\
                Number("rate",
                       EVENT_NOEVENT, xpos+BUTTON_WIDTH_AND_DX,\
                       ypos, NUMBER_WIDTH, NUMBER_HEIGHT,\
                       gparams["coarse_dense"]["rate"].val, 0.001, 4.0,
                       "The rate for coarsening dense areas")
    gparams["coarse_dense"]["numiter"] = \
                Number("Num iter",
                       EVENT_NOEVENT,
                       xpos+BUTTON_WIDTH_AND_DX+NUMBER_WIDTH_AND_DX,
                       ypos, NUMBER_WIDTH, NUMBER_HEIGHT,
                       gparams["coarse_dense"]["numiter"].val, 1, 15,
                       "The number of iteration for coarsening dense areas")

def draw_coarse_flat(xpos, ypos):
    global gparams
    Button("Coarse Flat", EVENT_COARSE_FLAT, xpos, ypos,
           BUTTON_WIDTH, BUTTON_HEIGHT, "Coarse flat areas of the mesh")
    gparams["coarse_flat"]["rate"] = \
            Number("rate", EVENT_NOEVENT,
                   xpos+BUTTON_WIDTH_AND_DX,
                   ypos, NUMBER_WIDTH, NUMBER_HEIGHT,
                   gparams["coarse_flat"]["rate"].val, 0.00001, 0.5,
                   "The rate for coarsening flat areas")

def update_min_max_angle():
    global gparams
    gparams["smooth"]["min_max_angle"].val = 180 - \
                2*gparams["smooth"]["max_min_angle"].val

def draw_smooth(xpos, ypos):
    global gparams
    Button("Smooth", EVENT_SMOOTH, xpos, ypos,
           BUTTON_WIDTH, BUTTON_HEIGHT, "Smooth the mesh")
    gparams["smooth"]["max_min_angle"] = \
            Number("Max min angle",
                   EVENT_UPDATE_MIN_MAX_ANGLE, xpos+BUTTON_WIDTH_AND_DX,
                   ypos, NUMBER_WIDTH, NUMBER_HEIGHT,
                   gparams["smooth"]["max_min_angle"].val, 10, 20,
                   "The maximal minimum angle for smoothing")
    gparams["smooth"]["max_iter"] = \
            Number("Max iter", EVENT_NOEVENT,
                   xpos+BUTTON_WIDTH_AND_DX+NUMBER_WIDTH_AND_DX,
                   ypos, NUMBER_WIDTH, NUMBER_HEIGHT,
                   gparams["smooth"]["max_iter"].val, 1, 30,
                   "The number of iteration for smooting mesh")
    Toggle("preserve ridges", EVENT_TOGGLE_PRESERVE_RIDGES,
           xpos+BUTTON_WIDTH_AND_DX+NUMBER_WIDTH_AND_DX*2,
           ypos, BUTTON_WIDTH, BUTTON_HEIGHT, gparams["smooth"]["preserve_ridges"].val,
           "preserve ridges during smoothing")

def draw_materials(xpos, ypos):
    global markers, material_button, marker_value
    # FIXME: Need to get the materials from the global scene!
    material_names, materials = get_materials()
    if materials:
        for mtrl in materials:
            if mtrl not in markers:
                markers[mtrl] = MARKER_DEFAULT
    if len(materials) != len(markers):
        for mtrl in markers.keys():
            if mtrl not in materials:
                markers.pop(mtrl)
    names = "|".join("%s %%x%d"%(name, i)
             for i, name in enumerate(material_names))
    material_button = Menu(names, EVENT_MATERIAL_SET, xpos, ypos,
                   BUTTON_WIDTH, BUTTON_HEIGHT, material_button.val,
                   "Choose a material.")
    marker_value = Number("Marker",
        EVENT_UPDATE_MARKER_VALUE, xpos+BUTTON_WIDTH_AND_DX,
        ypos, NUMBER_WIDTH, NUMBER_HEIGHT,
        marker_value.val, 1, 3000,
        "The marker value assigned to a boundary")

def draw_tetrahedralize(xpos, ypos):
    global tetmesh_value, tetmesh_name
    Button("Tetrahedralize", EVENT_TETRAHEDRALIZE, xpos, ypos,
           BUTTON_WIDTH, BUTTON_HEIGHT, "Tetrahedralize a surface mesh")
    
    #tetmesh_name = String("", EVENT_NOEVENT,
    #              xpos+BUTTON_WIDTH_AND_DX+NUMBER_WIDTH_AND_DX,
    #              ypos, NUMBER_WIDTH, NUMBER_HEIGHT,
    #              tetmesh_name.val,
    #              30, "Name of the tetrahedralized mesh")

Register(draw, event, bevent)

######################################################
# Main Body
######################################################
def toggle_overwrite():
    global overwrite_mesh
    overwrite_mesh = not overwrite_mesh

def toggle_preserve_ridges():
    global gparams
    gparams["smooth"]["preserve_ridges"].val = 0 if gparams["smooth"]["preserve_ridges"].val else 1

def repair_mesh():
    scn = Blender.Scene.GetCurrent()
    obj = scn.objects.active
    if not obj or obj.type != 'Mesh':
        BPyMessages.Error_NoMeshActive()
        return None, None
    
    # Get Blender mesh
    Blender.Window.WaitCursor(1)
    bmesh = obj.getData(mesh=1)

    # Collect quads
    nquads = 0

    # Remove free vertices
    vert_users = np.zeros(len(bmesh.verts))
    for f in bmesh.faces:
        f.sel = 0
        for v in f:
            vert_users[v.index] += 1
        if len(f) == 4:
            nquads += 1
            f.sel = 1
            
    if nquads:
        print "Found %d quads"%nquads

    for e in bmesh.edges:
        for v in e: 
            vert_users[v.index] += 1
    
    verts_free = (vert_users==0).nonzero()[0].tolist()
	
    if verts_free:
        print "Removed %s vertices"%len(verts_free)
        bmesh.verts.delete(verts_free)
    
    # Remove edges with no face connected to it
    edges = set() 
	
    for f in bmesh.faces:
        for edkey in f.edge_keys:
            edges.add(edkey)
	
    edges_free = []
    for e in bmesh.edges:
        if e.key not in edges:
            edges_free.append(e)
    if edges_free:
        bmesh.edges.delete(edges_free)
        print "Removed %s edges"%len(edges_free)

    # Convert selected quads to triangles
    bmesh.quadToTriangle()

    # Check for non-manifolds and open edges
    edge_map = dict((edge.key, []) for edge in bmesh.edges)
    
    for face in bmesh.faces:
        # Unselect all faces
        face.sel = 0
        for edge in face.edge_keys:
            edge_map[tuple(sorted(edge))].append(face.index)
    
    open_edges = []
    non_manifold_edges = []
    open_vertices_map = {}
    for edge, faces in edge_map.iteritems():
        if len(faces) == 1:
            open_edges.append(edge)
            for i, vert_id in enumerate(edge):
                if vert_id not in open_vertices_map:
                    open_vertices_map[vert_id] = [edge[1-i]]
                else:
                    open_vertices_map[vert_id].append(edge[1-i])
                    open_vertices_map[vert_id].sort()
        if len(faces) > 2:
            non_manifold_edges.append(edge)
    
    print "Found %d open edges"%len(open_edges)
    print "Found %d non manifold edges"%len(non_manifold_edges)

    free_faces = []
    add_faces = []
    for edge in non_manifold_edges:
        for face in edge_map[edge]:
            open_face_edges = [face_edge for face_edge in \
                               bmesh.faces[face].edge_keys\
                               if face_edge in open_edges]
            if len(open_face_edges) == 2:
                print "Found a totally open face"
                free_faces.append(face)
                for open_edge in open_face_edges:
                    open_edges.remove(open_edge)
                    for vert_id in open_edge:
                        if len(open_vertices_map[vert_id]) == 1:
                            open_vertices_map.pop(vert_id)
                        else:
                            open_vertices_map[vert_id].pop\
                                        (1-open_edge.index(vert_id))
            if len(open_face_edges) == 1:
                bmesh.faces[face].sel = 1
                print "Found a complex connected face. Selects it!"

    #for vert, edges in open_vertices_map:
    #print open_vertices_map
    open_vertices = sorted(open_vertices_map.keys())
    for vert in open_vertices:
        if vert in open_vertices_map and \
           open_vertices_map[vert][0] in open_vertices_map and \
           open_vertices_map[vert][1] in open_vertices_map and \
           open_vertices_map[open_vertices_map[vert][0]][0] == vert and \
           open_vertices_map[open_vertices_map[vert][1]][0] == vert:
            add_faces.append([vert]+open_vertices_map[vert])
            for vert_ind in open_vertices_map[vert]:
                open_vertices_map.pop(vert_ind)
            open_vertices_map.pop(vert)

    # If not connected all faces
    if open_vertices_map:
        print "Found a non trivial connected open edges selects face"
        for vert0, verts in open_vertices_map.iteritems():
            edges = [tuple(sorted([vert0, vert1])) for vert1 in verts]
            for edge in edges:
                bmesh.faces[edge_map[edge][0]].sel = 1
    
    # Free and add faces
    bmesh.faces.delete(1, free_faces)
    bmesh.faces.extend(add_faces)
            
    # Remove free vertices
    vert_users = np.zeros(len(bmesh.verts))
    for f in bmesh.faces:
        for v in f:
            vert_users[v.index] += 1

    for e in bmesh.edges:
        for v in e: 
            vert_users[v.index] += 1
    
    verts_free = (vert_users==0).nonzero()[0].tolist()
	
    if verts_free:
        print "Removed %s vertices"%len(verts_free)
        bmesh.verts.delete(verts_free)

    # Harmonize the normals
    for face in bmesh.faces:
        face.sel = 1
    print "Recalculate normals"
    bmesh.recalcNormals(0)
    

# Populate the global name space with GAMer action functions
def blender_to_gamer():
    "Transfer the active mesh to a GAMer surface mesh"
    scn = Blender.Scene.GetCurrent()
    obj = scn.objects.active
    if not obj or obj.type != 'Mesh':
        BPyMessages.Error_NoMeshActive()
        return None, None
    
    # Get Blender mesh
    Blender.Window.WaitCursor(1)
    bmesh = obj.getData(mesh=1)

    # Transfere data from blender mesh to gamer mesh
    gmesh = gamer.SurfaceMesh(len(bmesh.verts),\
                              len(bmesh.faces))
    
    for i, bverts in enumerate(bmesh.verts):
        gvert = gmesh.vertex(i)
        gvert.x, gvert.y, gvert.z = tuple(bverts.co)

    for i, bface in enumerate(bmesh.faces):
        gface = gmesh.face(i)
        gface.a, gface.b, gface.c = \
                 bface.v[0].index, bface.v[1].index, bface.v[2].index
    
    return bmesh.name, gmesh

def gamer_to_blender(gmesh, mesh_name="gamer_improved"):
    # Check arguments
    if not isinstance(gmesh, gamer.SurfaceMesh):
        raise TypeError, "expected a SurfaceMesh"

    verts = []
    faces = []
    for i in xrange(gmesh.num_vertices):
        gvert = gmesh.vertex(i)
        verts.append((gvert.x, gvert.y, gvert.z))

    for i in xrange(gmesh.num_faces):
        gface = gmesh.face(i)
        faces.append((gface.a, gface.b, gface.c))

    # Get scene
    scn = Blender.Scene.GetCurrent()

    # FIXME: Does not work...
    if overwrite_mesh:
        
        obj = scn.objects.active
        if not obj or obj.type != 'Mesh':
            BPyMessages.Error_NoMeshActive()
            return 
        Blender.Window.WaitCursor(1)
        bmesh = obj.getData(mesh=1)

        # Update the coordinates
        for gv, bv in zip(verts, bmesh.verts):
            bv.co.x, bv.co.y, bv.co.z = gv

        # Update the faces
        for gf, bf in zip(faces, bmesh.faces):
            bf.verts = tuple([bmesh.verts[i] for i in gf])

        # Remove superflous faces
        if len(bmesh.faces) > len(faces):
            print range(len(faces), len(bmesh.faces))
            bmesh.faces.delete(1, range(len(faces), len(bmesh.faces)))
        # Extend additional faces
        elif len(bmesh.faces) < len(faces):
            bmesh.faces.extend(faces[range(len(faces)-1, len(bmesh.faces)-1, -1)])

        # Remove superflous vertices
        if len(bmesh.verts) > len(verts):
            print range(len(verts), len(bmesh.verts))
            bmesh.verts.delete(range(len(verts), len(bmesh.verts)))
        # Extend additional vertices
        elif len(bmesh.verts) < len(verts):
            bmesh.verts.extend(verts[range(len(verts)-1, len(bmesh.verts)-1, -1)])

    else:
        # Switch to another layer
        switch_to_layer = scn.getLayers()[-1] + 1
        switch_to_layer = 1 if switch_to_layer > 20 else switch_to_layer
        scn.setLayers([switch_to_layer])

        # Create new mesh
        bmesh = Blender.Mesh.New(mesh_name)
        bmesh.verts.extend(verts)
        bmesh.faces.extend(faces)

        # Apply 
        scn.objects.active = scn.objects.new(bmesh, mesh_name)
    Blender.Window.RedrawAll()

def tetrahedralize_action(filename):
    "Callback function for the tetrahedralize File chooser"

    # Load options and materials from registry
    load_from_registry()
    
    # Get gamer mesh
    name, gmesh = blender_to_gamer()
    
    # Get materials on mesh
    names, mtrls = get_materials()
    
    # If there are materials defined in the mesh convert the marker information
    # to the gamer mesh
    if mtrls:
        gmesh.reset_face_markers()
        scn = Blender.Scene.GetCurrent()
        obj = scn.objects.active
        if not obj or obj.type != 'Mesh':
            return

        Blender.Window.WaitCursor(1)
        bmesh = obj.getData(mesh=1)
        for i, bf in enumerate(bmesh.faces):
            gmesh.set_face_marker(i, markers[bmesh.materials[bf.mat]])

    # Tetrahedralize mesh
    gem_mesh = gamer.GemMesh(gmesh)
    
    # Write to file
    if tetmesh_value.val == 0:
        gem_mesh.write_mcsf(filename)

def tetrahedralize_mesh():
    global tetmesh_value

    names = "Choose output format%t|" + "|".join("%s %%x%d"%(name, i)
                            for i, name in enumerate(TETMESH_OPTIONS))
    tetmesh_value.val = PupMenu(names) 

    if tetmesh_value.val < 0:
        return
    
    # Save any options and materials to registry
    save_to_registry()

    # Call the File selector
    Blender.Window.FileSelector(tetrahedralize_action, 'Tetrahedralize',\
                                "*"+TETMESH_SUFFIX[tetmesh_value.val])

def print_eigenvalues():
    name, gmesh = blender_to_gamer()
    if gmesh is None:
        return
    gmesh.eigenvalues()
    
def gamer_action(action):
    "A function that returns gamer action functions"
    def action_func():
        kwargs = dict()
        for name, value in gparams[action].iteritems():
            kwargs[name] = value.val
        name, gmesh = blender_to_gamer()
        if gmesh is None:
            return
        if GAMER_SUFFIX in name:
            base_name = name.split(GAMER_SUFFIX)[0]
            next_tick = int(re.findall(GAMER_SUFFIX+"([0-9].)", name)[0])+1
            new_name = base_name + GAMER_SUFFIX + ("%%0%dd"%NUM_TICK_FIGURES)%next_tick
        else:
            new_name_length = len(name) + len(GAMER_SUFFIX) + NUM_TICK_FIGURES
            if new_name_length>MAX_MESH_LENGTH_NAME:
                name = name[:MAX_MESH_LENGTH_NAME-new_name_length]
            new_name = name + GAMER_SUFFIX + "0"*NUM_TICK_FIGURES
        if gmesh is None:
            return
        Redraw(1)
        getattr(gmesh, action)(**kwargs)
        gamer_to_blender(gmesh, new_name)
    return action_func

# Fill the event dict with callback functions
event_dict = {EVENT_EXIT:exit_event,
              EVENT_UPDATE_MIN_MAX_ANGLE:update_min_max_angle,
              EVENT_TOGGLE_OVERWRITE:toggle_overwrite,
              EVENT_TOGGLE_PRESERVE_RIDGES:toggle_preserve_ridges,
              EVENT_MATERIAL_SET:select_material,
              EVENT_UPDATE_MARKER_VALUE:update_marker_value,
              EVENT_TETRAHEDRALIZE:tetrahedralize_mesh,
              EVENT_EIGENVALUES:print_eigenvalues,
              EVENT_REPAIR_MESH:repair_mesh,
              }

for action in gparams:
    try:
        exec("event_dict[%s] = gamer_action('%s')"%("EVENT_%s"%(action.upper()), action))
    except:
        pass

