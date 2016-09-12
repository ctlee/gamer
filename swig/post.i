/* -*- C -*- */
// Copyright (C) 2010 Johan Hake
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2010-08-06
// Last changed: 2013-03-06
//
// ===========================================================================
// SWIG directives for GAMer 
//
// The directives in this file are applied _after_ the header files has been 
// loaded.
// ===========================================================================

// Code that will both be verbitim inserted in the header section in the 
// generated wrapper code, but also parsed by SWIG to generate wrapper code
%inline%{
    
static FLTVECT* _Vertex_getitem(SurfaceMesh* surfmesh, unsigned int i){
  return surfmesh->vertex + i;
}

static INT3VECT* _Face_getitem(SurfaceMesh* surfmesh, unsigned int i){
  return surfmesh->face + i;
}

static FETK_VX* _GemMesh_vertex_getitem(GemMesh* gem_mesh, unsigned int i){
  return gem_mesh->vv + i;
}

static FETK_SS* _GemMesh_tet_getitem(GemMesh* gem_mesh, unsigned int i){
  return gem_mesh->ss + i;
}

%}

// Extend both the C and Python interface of SurfaceMesh 
%extend SurfaceMesh {

  // Create a surface mesh by triangulating a list of connected segments
  SurfaceMesh(double* pointlist, int numberofcoordinates, 
	      char* triangle_quality_params="q20")
  {
    char triangle_params[50];
    char default_params[] = "pYz"; //Y
    sprintf(triangle_params, "%s%s", default_params, triangle_quality_params);
    return SurfaceMesh_triangulate(pointlist, numberofcoordinates, triangle_params);
  }

  // Create an empty mesh 
  SurfaceMesh(unsigned int num_vertices, unsigned int num_faces)
  {
    return SurfaceMesh_ctor(num_vertices, num_faces);
  }
  
  // Create a spherical surface mesh
  SurfaceMesh(unsigned int level)
  {
    return SurfaceMesh_sphere(level);
  }
  
  // Constructor for pdb or off files
  SurfaceMesh(char* filename)
  {
    char* basename_end_pointer;
  
    // Check if input is a .off file
    basename_end_pointer = strstr(filename, ".off");
    if (basename_end_pointer != NULL){
      return SurfaceMesh_readOFF(filename);
    }
  
    // Check if input is a .poly file
    basename_end_pointer = strstr(filename, ".poly");
    if (basename_end_pointer != NULL)
    {
      char basename[100];
      basename[0] = '\0';

      // Pointer arithmetic
      strncat(basename, filename, basename_end_pointer-filename);
      return SurfaceMesh_readPoly(basename);
    }
  
    // If not correct file format
    printf("Provide a \'.off\' or \'.poly\'file.\n");
    return SurfaceMesh_ctor(0,0);
  }
  
  // Constructor for LatticeData
  SurfaceMesh(char* segmentation_filename, float isovalue, bool padd)
  {
    return SurfaceMesh_readLattice(segmentation_filename, isovalue, padd);
  }  
  
  // Constructor for LatticeData together with intensity values
  SurfaceMesh(char* segmentation_filename, char* intensity_filename, 
              float isovalue, bool padd)
  {
    return SurfaceMesh_readLattice(segmentation_filename, 
        			   intensity_filename, isovalue, padd);
  }  
  
  // Constructor that merges two surface meshes
  SurfaceMesh(SurfaceMesh* in0, SurfaceMesh* in1)
  {
    return SurfaceMesh_merge(in0, in1);
  }
  
  // Destructor
  ~SurfaceMesh()
  {
    SurfaceMesh_dtor($self);
  }
  
  // Python interface
  %pythoncode %{
    def vertex(self, i):
        "Return the ith vertex"
        if not isinstance(i, int):
            raise TypeError("expected an int for the index value")
        if i >= self.num_vertices:
            raise IndexError("index out of range")
        return _Vertex_getitem(self, i)
    
    def vertices(self):
        "Return an iterator over vertices"
        for i in range(self.num_vertices):
            yield _Vertex_getitem(self, i)
    
    def face(self, i):
        "Return the ith vertex"
        if not isinstance(i, int):
            raise TypeError("expected an int for the index value")
        if i >= self.num_faces:
            raise IndexError("index out of range")
        return _Face_getitem(self, i)
    
    def faces(self):
        "Return an iterator over faces"
        for i in range(self.num_faces):
            yield _Face_getitem(self, i)
    
    def refine(self):
        "Refine a surface mesh uniformly"
        SurfaceMesh_refine(self)
    
    def smooth(self, max_min_angle=15, min_max_angle=150,
               max_iter=6, preserve_ridges=False):
        "Smooth the mesh"
        return SurfaceMesh_smooth(self, max_min_angle, min_max_angle, 
        			  max_iter, preserve_ridges)
    
    def normal_smooth(self):
        "Smooth the facet normals of the mesh"
        SurfaceMesh_normalSmooth(self)
    
    def correct_normals(self):
        "Correct any normals that might be screwed"
        SurfaceMesh_correctNormals(self)
    
    def flip_normals(self):
        "Flip the direction of the normals"
        SurfaceMesh_flipNormals(self)
    
    def write_off(self, filename):
        "Write the mesh to file"
        SurfaceMesh_writeOFF(self, filename)
    
    def write_poly(self, filename):
        "Write the mesh to a poly file"
        SurfaceMesh_writePoly(self, filename)
    
    def coarse_flat(self, rate=0.016):
        "Coarse flat areas"
        SurfaceMesh_coarse(self, rate, 1, 0, -1)

    def scale(self, *args):
        """
        Scale all coordiantes

        Arguments:
           scale : Single float
           x_scale, y_scale, z_scale : a value for each coordinate
        """
        SurfaceMesh_scale(self, *args)

    def centralize(self):
        "Centralize the mesh"
        SurfaceMesh_centeralize(self)

    def translate(self, dx, dy, dz):
        "Translate mesh"
        SurfaceMesh_translate(self, dx, dy, dz)
        
    def coarse_dense(self, rate=1.6, numiter=1):
        "Coarse flat areas"
        if not isinstance(rate, (int, float)) or rate <= 0:
            raise TypeError("expected a positive scalar for the 'rate' argument")
        if not isinstance(numiter, int) or numiter < 1:
            raise TypeError("expected a positive scalar for the 'numiter' argument")
        for i in range(numiter):
            SurfaceMesh_coarse(self, rate, 0, 10, -1)
    
    def get_center_radius(self):
        "Return an ATOM structor of the radius and center of a surface mesh"
        return SurfaceMesh_getCenterRadius(self)

    def remove_unconnected_patches(self, minimal_number):
        "Remove patches of faces with size smaller than minimal_number"
        SurfaceMesh_removeUnconnectedPatches(self, minimal_number)
  %}
  
  //    def getMinMaxAngles(self, max_min_angle=15, min_max_angle=150):
  //        """Return angle statistics"""
  //        return SurfaceMesh_getMinMaxAngles(max_min_angle, min_max_angle)
  
}

%pythoncode%{
def read_PDB_gauss(filename, blobbyness=-0.2, iso_value=2.5):
    "Read a PDB file and generate a SurfaceMesh using Gaussian iso surfaces"
    return SurfaceMesh_readPDB_gauss(filename, blobbyness, iso_value)
%}


// Extend both the C and Python interface of GemMesh
%extend GemMesh {

  GemMesh(unsigned int num_vertices, unsigned int num_tetrahedrons, bool higher_order=false)
  {
    return GemMesh_ctor(num_vertices, num_tetrahedrons, higher_order);
  }
  
  GemMesh(SurfaceMesh** surfmeshes, int num_meshes, char* tetgen_quality_params="q2.0qq18")
  {
    char tetgen_params[50];
    char default_params[] = "KnnepAAV"; //Y
    sprintf(tetgen_params, "%s%s", default_params, tetgen_quality_params);
    return GemMesh_fromSurfaceMeshes(surfmeshes, num_meshes, tetgen_params);
  }
  
  GemMesh(SurfaceMesh* surfmesh, char* tetgen_quality_params="q2.0qq18")
  {
    SurfaceMesh* surfmeshes[1] = {surfmesh};
    return GemMesh_fromSurfaceMeshes(surfmeshes, 1, tetgen_quality_params);
  }
  
  ~GemMesh()
  {
    GemMesh_dtor($self);
  }
  
  void write_mcsf(char* filename)
  {
    GemMesh_writeMcsf(self, filename);
  }
  
  void write_dolfin(char* filename)
  {
    GemMesh_writeDolfin(self, filename);
  }

  void write_carp(char* filename, PyObject *op)
  {
    if (!PyList_Check(op))
      throw std::runtime_error("Expected a list of strings as second argument.");
    int i = 0;
    int argc = PyList_Size(op);
    char **argv = (char **) malloc((argc+1)*sizeof(char *));
    for (i = 0; i < argc; i++)
    {
      PyObject *o = PyList_GetItem(op,i);
      if (PyString_Check(o))
	argv[i] = PyString_AsString(o);
      else
      {
	free(argv);
	throw std::runtime_error("Expected a list of strings as second argument.");
      }
    }
    argv[i] = 0;
    GemMesh_writeCarp(self, filename, argc, argv);
    free(argv);
  }
  
  void write_diffpack(char* filename, PyObject *op)
  {
    if (!PyList_Check(op))
      throw std::runtime_error("Expected a list of strings as second argument.");
    int i = 0;
    int argc = PyList_Size(op);
    char **argv = (char **) malloc((argc+1)*sizeof(char *));
    for (i = 0; i < argc; i++)
    {
      PyObject *o = PyList_GetItem(op,i);
      if (PyString_Check(o))
	argv[i] = PyString_AsString(o);
      else
      {
	free(argv);
	throw std::runtime_error("Expected a list of strings as second argument.");
      }
    }
    argv[i] = 0;
    GemMesh_writeDiffpack(self, filename, argc, argv);
    free(argv);
  }
  
  void write_off(char* filename)
  {
    GemMesh_writeOFF(self, filename);
  }
  
  %pythoncode%{
    def vertex(self, i):
        "Return the ith vertex"
        if not isinstance(i, int):
            raise TypeError("expected an int for the index value")
        if i >= self.num_vertices:
            raise IndexError("index out of range")
        return _GemMesh_vertex_getitem(self, i)
    
    def vertices(self):
        "Return an iterator over vertices"
        for i in range(self.num_vertices):
            yield _GemMesh_vertex_getitem(self, i)
    
    def cell(self, i):
        "Return the ith cell"
        if not isinstance(i, int):
            raise TypeError("expected an int for the index value")
        if i >= self.num_cells:
            raise IndexError("index out of range")
        return _GemMesh_tet_getitem(self, i)
    
    def cells(self):
        "Return an iterator over faces"
        for i in range(self.num_cells):
            yield _GemMesh_tet_getitem(self, i)
	      
    %}

  }
