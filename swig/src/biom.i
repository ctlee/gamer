/**
 *  @file       biom.h
 *  @ingroup    global_gamer
 *  @brief      Some parameter/datatype definitions
 *  @author     Zeyun Yu (zeyun.yu@gmail.com)
 *  @note       None
 *  @version    $Id: biom.i,v 1.8 2013/03/06 22:55:13 fetk Exp $
 *  
 *  @attention
 *  @verbatim
 *
 * GAMER = < Geometry-preserving Adaptive MeshER >
 * Copyright (C) 1994-- Michael Holst and Zeyun Yu
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 *  @endverbatim
 */

// critical parameters for mesh generation

/** @brief Isovalue used in the Marching Cube method */
#define IsoValue          2.5      
/** @brief Blurring blobyness used in conversion from PDB/PQR to 3D volumes */
#define BLOBBYNESS        -0.2f    
/** @brief Discretization rate of 3D volumes */
#define DIM_SCALE         1.99   
/** @brief Coarsening Rate in surface post-processing */
#define CoarsenRate       0.1666   
/** @brief The minimal volumes (in voxels) of islands to be removed */
#define MIN_VOLUME        333333   
/** @brief The size of the bounding sphere (= object size X the following rate) */
#define SphereRatio       40

#define MaxVal               999999
#define MaxAtom              10

// Other definitions and data structures

/** @brief Other definition */
#define PIE              3.14159265358979f 
/** @brief Other definition */
#define IndexVect(i,j,k) ((k)*xdim*ydim + (j)*xdim + (i))
#define IndexVect1(i,j,k) ((k)*xdim1*ydim1 + (j)*xdim1 + (i))
/** @brief Other definition */
#define max(x, y)        ((x>y) ? (x):(y))
/** @brief Other definition */
#define min(x, y)        ((x<y) ? (x):(y))

/** @brief Other definition */
#define _LITTLE_ENDIAN   1

/** @brief Other data structure FLT2VECT (float) */
typedef struct 
{
  float x;   /**< @brief x-coordinate */
  float y;   /**< @brief y-coordinate */
}FLT2VECT;

/** @brief Other data structure FLTVECT (float) */
typedef struct 
{
  float x;   /**< @brief x-coordinate */
  float y;   /**< @brief y-coordinate */
  float z;   /**< @brief z-coordinate */
  int m;     /**< @brief Marker */
  bool sel;  /**< @brief selection flag */
}FLTVECT;

/** @brief Other data structure DBLVECT (double float) */
typedef struct 
{
  double x;   /**< @brief x-coordinate */
  double y;   /**< @brief y-coordinate */
  double z;   /**< @brief z-coordinate */
}DBLVECT;

/** @brief Other data structure INT3VECT (int) */
typedef struct 
{
  int a;     /**< @brief first integer */
  int b;     /**< @brief second integer */
  int c;     /**< @brief third integer */
  int m;     /**< @brief Marker */
  bool sel;  /**< @brief selection flag */
}INT3VECT;

/** @brief Other data structure INT4VECT (int) */
typedef struct 
{
  int a;   /**< @brief first integer */
  int b;   /**< @brief second integer */
  int c;   /**< @brief third integer */
  int d;   /**< @brief fourth integer */
}INT4VECT;

/** @brief Other data structure NPNT3 */
typedef struct NeighborPoint3 NPNT3;
/** @brief Other data structure NeighborPoint3 */
struct NeighborPoint3
{
  int a;   /**< @brief first integer */
  int b;   /**< @brief second integer */
  int c;   /**< @brief third integer */
  NPNT3 *next;   /**< @brief pointer to the next triangle */
};

/** @brief Other data structure SurfaceMesh (for surface meshes) */
typedef struct 
{
  int num_vertices;   /**< @brief number of vertices */
  int num_faces;      /**< @brief number of triangles */
  int marker;            /**< @brief doman marker, to be used when tetrahedralizing */
  float volume_constraint; /**< @brief volume constraint of the tetrahedralized domain */
  bool use_volume_constraint; /**< @brief flag that determines if the volume constraint is used */
  bool as_hole; /**< @brief flag that determines if the mesh is a hole or not */
}SurfaceMesh;

typedef struct FETK_VX
{
  int id;   
  int chrt; 
  float x;  
  float y;  
  float z;  
}FETK_VX;

typedef struct FETK_SS
{
  int id;  
  int grp;
  int mat;
  int fa;  
  int fb;  
  int fc;  
  int fd;  
  int na;  
  int nb;  
  int nc;  
  int nd;  
}FETK_SS;

typedef struct GemMesh
{
  int dim;       
  int dimii;     
  int num_vertices;
  int num_cells;
  bool higher_order;
}GemMesh;


/** @brief Other data structure SPNT */
typedef struct SamplePoint SPNT;
/** @brief Other data structure SamplePoint */
struct SamplePoint
{
  float x;   /**< @brief x-coordinate */
  float y;   /**< @brief y-coordinate */
  float z;   /**< @brief z-coordinate */
  SPNT *next;   /**< @brief pointer to next vertex */
};


/** @brief Other data structure NPNT2 */
typedef struct NeighborPoint2 NPNT2;
/** @brief Other data structure NeighborPoint2 */
struct NeighborPoint2
{
  int a;   /**< @brief first integer */
  int b;   /**< @brief second integer */
  NPNT2 *next;   /**< @brief pointer to the next point */
};

/** @brief Other data structure ATOM */
typedef struct 
{
  float x;   /**< @brief x-coordinate */
  float y;   /**< @brief y-coordinate */
  float z;   /**< @brief z-coordinate */
  float radius;   /**< @brief radius */
}ATOM;


/** @brief Other data structure EIGENVECT */
typedef struct 
{
  float x1;   /**< @brief x-coordinate of first eigenvector */
  float y1;   /**< @brief y-coordinate of first eigenvector */
  float z1;   /**< @brief z-coordinate of first eigenvector */
  float x2;   /**< @brief x-coordinate of second eigenvector */
  float y2;   /**< @brief y-coordinate of second eigenvector */
  float z2;   /**< @brief z-coordinate of second eigenvector */
  float x3;   /**< @brief x-coordinate of third eigenvector */
  float y3;   /**< @brief y-coordinate of third eigenvector */
  float z3;   /**< @brief z-coordinate of third eigenvector */
}EIGENVECT;


/** @brief Other data structure MinHeapS */
typedef struct 
{
  unsigned short *x; /**< @brief x-coordinate */
  unsigned short *y; /**< @brief y-coordinate */
  unsigned short *z; /**< @brief z-coordinate */
  int *seed;         /**< @brief seed */
  float *dist;       /**< @brief distance */
  int size;          /**< @brief size */
}MinHeapS;


/** @brief Other data structure SEEDS */
typedef struct
{
  float seedx;       /**< @brief x-coordinate */
  float seedy;       /**< @brief y-coordinate */
  float seedz;       /**< @brief z-coordinate */
  int atom[MaxAtom]; /**< @brief atom array */
}SEEDS;

// GemMesh functions
//GemMesh* GemMesh_fromTetgen(tetgenio& tetio);
//GemMesh* GemMesh_fromSurfaceMesh(SurfaceMesh* surfmesh, char* tetgen_params);
//GemMesh* GemMesh_fromPdb(tetgenio* out, float radius, float centerx, 
//			 float centery, float centerz, char *ActiveSite, 
//			 int output_flag);
void GemMesh_writeMcsf(GemMesh* out, char* filename);
void GemMesh_dtor(GemMesh*);
void GemMesh_writeOFF(GemMesh* Gem_mesh, char* filename);

// Not working...
//void SurfaceMesh_getMinMaxAngles(SurfaceMesh* ,float*, float*, int*, int*, int, int);
void SurfaceMesh_refine(SurfaceMesh*);
bool SurfaceMesh_smooth(SurfaceMesh* surfmesh, 
			unsigned int max_min_angle, 
			unsigned int min_max_angle, 
			unsigned int max_iter,
			bool flip_edges);
void SurfaceMesh_normalSmooth(SurfaceMesh* surfmesh);
char SurfaceMesh_coarse(SurfaceMesh* surfmesh, float coarse_rate,
			float flatness_rate, float denseness_weight,
			float max_normal_angle);
ATOM SurfaceMesh_getCenterRadius(SurfaceMesh* surfmesh);
// Not working
// void SurfaceMesh_assignActiveSites(SurfaceMesh*, ATOM*, unsigned int, unsigned int*, unsigned int*);
ATOM SurfaceMesh_getCenterRadius(SurfaceMesh* surfmesh);
void SurfaceMesh_translate(SurfaceMesh* surfmesh, float dx, float dy, float dz);
void SurfaceMesh_scale(SurfaceMesh* surfmesh, float scale_x, 
		       float scale_y, float scale_z);
void SurfaceMesh_scale(SurfaceMesh* surfmesh, float scale);
void SurfaceMesh_centeralize(SurfaceMesh* surfmesh);
void SurfaceMesh_flipNormals(SurfaceMesh* surfmesh);
void SurfaceMesh_correctNormals(SurfaceMesh* surfmesh);
void SurfaceMesh_writeOFF(SurfaceMesh* surfmesh, char* filename);
void SurfaceMesh_writePoly(SurfaceMesh* surfmesh, char* filename);
void SurfaceMesh_removeUnconnectedPatches(SurfaceMesh* surfmesh, int minimal_number);
void SurfaceMesh_removeUnconnectedVertices(SurfaceMesh* surfmesh);
void SurfaceMesh_deleteVertices(SurfaceMesh* surfmesh);
void SurfaceMesh_deleteFaces(SurfaceMesh* surfmesh);

SurfaceMesh* SurfaceMesh_readPDB_molsurf(char* filename);
SurfaceMesh* SurfaceMesh_readPDB_gauss(char* filename, float blobbyness, 
				       float iso_value);

// void ReadActiveSiteFile(char*, unsigned int&, ATOM*&, unsigned int*&, float*&);
void ReadActiveSiteFile(char*, unsigned int&, ATOM*&, unsigned int*&);
void ReadRawiv(int *, int *, int *, float **, char *,float *, float *);
void SurfaceExtract(TeTraMesh*, SurfaceMesh*);
// Not working
// float PDB2Volume(char *, float **, int *, int *, int *,float *, float *, ATOM **,int *,char);
// void MarchingCube(int, int, int, float *, float, float, SurfaceMesh **,SPNT **);
