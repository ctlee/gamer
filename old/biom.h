/**
 *  @file       biom.h
 *  @ingroup    global_gamer
 *  @brief      Some parameter/datatype definitions
 *  @author     Zeyun Yu (zeyun.yu@gmail.com)
 *  @note       None
 *  @version    $Id: biom.h,v 1.38 2013/03/06 22:55:13 fetk Exp $
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

#ifndef _BIOM_H_
#define _BIOM_H_

#include <vector>
#include <tuple>

#include "gamer_base.h"
#include <libraries/triangle/triangle.h>
#include <libraries/tetgen/tetgen.h>

#if 0
# include <stdio.h>     /* One of 15 ISO-C headers -- get via MALOC */
# include <stdlib.h>    /* One of 15 ISO-C headers -- get via MALOC */
# include <math.h>      /* One of 15 ISO-C headers -- get via MALOC */
# include <time.h>      /* One of 15 ISO-C headers -- get via MALOC */
# include <string.h>    /* One of 15 ISO-C headers -- get via MALOC */
# include <sys/times.h> /* System-dependent; include in .c file via MALOC */
# include <sys/stat.h>  /* System-dependent; include in .c file via MALOC */
# include <memory.h>    /* System-dependent; include in .c file via MALOC */
#endif // if 0

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
#define IndexVect(i, j, k) ((k) * xdim * ydim + (j) * xdim + (i))
#define IndexVect1(i, j, k) ((k) * xdim1 * ydim1 + (j) * xdim1 + (i))

/** @brief Other definition */

// #define max(x, y)        ((x>y) ? (x):(y))

/** @brief Other definition */

// #define min(x, y)        ((x<y) ? (x):(y))

/** @brief Other definition */
#define _LITTLE_ENDIAN   1

/** @brief Other data structure ATOM */
struct ATOM {
    float x;      /**< @brief x-coordinate */
    float y;      /**< @brief y-coordinate */
    float z;      /**< @brief z-coordinate */
    float radius; /**< @brief radius */
};

/** @brief Other data structure FLT2VECT (float) */
struct FLT2VECT {
    float x; /**< @brief x-coordinate */
    float y; /**< @brief y-coordinate */
};

/** @brief Other data structure FLTVECT (float) */
struct FLTVECT {
    float x;   /**< @brief x-coordinate */
    float y;   /**< @brief y-coordinate */
    float z;   /**< @brief z-coordinate */
    int   m;   /**< @brief Marker */
    bool  sel; /**< @brief selection flag */
};

/** @brief Other data structure INT3VECT (int) */
struct INT3VECT {
    int  a;   /**< @brief first integer */
    int  b;   /**< @brief second integer */
    int  c;   /**< @brief third integer */
    int  m;   /**< @brief Marker */
    bool sel; /**< @brief selection flag */
};

/** @brief Other data structure INT4VECT (int) */
struct INT4VECT {
    int a; /**< @brief first integer */
    int b; /**< @brief second integer */
    int c; /**< @brief third integer */
    int d; /**< @brief fourth integer */
};

/** @brief Other data structure INT6VECT (int) */
struct INT6VECT {
    int a; /**< @brief first integer */
    int b; /**< @brief second integer */
    int c; /**< @brief third integer */
    int d; /**< @brief fourth integer */
    int e; /**< @brief fifth integer */
    int f; /**< @brief sixth integer */
};

/** @brief Other data structure NPNT3 */
typedef struct NeighborPoint3 NPNT3;

/** @brief Other data structure NeighborPoint3 */
struct NeighborPoint3 {
    int    a;    /**< @brief first integer */
    int    b;    /**< @brief second integer */
    int    c;    /**< @brief third integer */
    NPNT3 *next; /**< @brief pointer to the next triangle */
};

/** @brief Other data structure Triface */
struct Triface {
    int      points[3]; /**< @brief point list */
    int      marker;    /**< @brief marker integer */
    Triface *next;      /**< @brief pointer to the next Triface */
};

/** @brief Other data structure TeTraMesh (for tetrahedral meshes) */
struct TeTraMesh {
    int       num_vertices;
    int       num_cells;
    FLTVECT  *vertex;   /**< @brief pointer to the vertices */
    INT4VECT *face;     /**< @brief pointer to the tetrahedra */
    INT4VECT *neighbor; /**< @brief pointer to the neighbors (tetrahedra) */
};

/** @brief Other data structure SPNT */
typedef struct SamplePoint SPNT;

/** @brief Other data structure SamplePoint */
struct SamplePoint {
    float x;    /**< @brief x-coordinate */
    float y;    /**< @brief y-coordinate */
    float z;    /**< @brief z-coordinate */
    SPNT *next; /**< @brief pointer to next vertex */
};

/** @brief Other data structure NPNT2 */
typedef struct NeighborPoint2 NPNT2;

/** @brief Other data structure NeighborPoint2 */
struct NeighborPoint2 {
    int    a;    /**< @brief first integer */
    int    b;    /**< @brief second integer */
    NPNT2 *next; /**< @brief pointer to the next point */
};


/** @brief Other data structure EIGENVECT */
struct EIGENVECT {
    float x1; /**< @brief x-coordinate of first eigenvector */
    float y1; /**< @brief y-coordinate of first eigenvector */
    float z1; /**< @brief z-coordinate of first eigenvector */
    float x2; /**< @brief x-coordinate of second eigenvector */
    float y2; /**< @brief y-coordinate of second eigenvector */
    float z2; /**< @brief z-coordinate of second eigenvector */
    float x3; /**< @brief x-coordinate of third eigenvector */
    float y3; /**< @brief y-coordinate of third eigenvector */
    float z3; /**< @brief z-coordinate of third eigenvector */
};


/** @brief Other data structure MinHeapS */
struct MinHeapS {
    unsigned short *x;    /**< @brief x-coordinate */
    unsigned short *y;    /**< @brief y-coordinate */
    unsigned short *z;    /**< @brief z-coordinate */
    int            *seed; /**< @brief seed */
    float          *dist; /**< @brief distance */
    int             size; /**< @brief size */
};


/** @brief Other data structure SEEDS */
struct SEEDS {
    float seedx;         /**< @brief x-coordinate */
    float seedy;         /**< @brief y-coordinate */
    float seedz;         /**< @brief z-coordinate */
    int   atom[MaxAtom]; /**< @brief atom array */
};

/** @brief Other data structure FET_VX (Vertex coordinates for tetrahedral mesh) */
struct FETK_VX {
    int   id;   /**< @brief id of the vertex */
    int   chrt; /**< @brief vertex marker */
    float x;    /**< @brief x coordinate of the vertex */
    float y;    /**< @brief y coordinate of the vertex */
    float z;    /**< @brief z coordinate of the vertex */
};

/** @brief Other data structure FET_SS (Cell data for tetrahedral mesh) */
struct FETK_SS {
    int id;  /**< @brief id of the cell */
    int grp; /**< @brief grp id of the cell */
    int mat; /**< @brief material marker */
    int fa;  /**< @brief First vertex id  */
    int fb;  /**< @brief Second vertex id  */
    int fc;  /**< @brief Third vertex id  */
    int fd;  /**< @brief Fourth vertex id  */
    int na;  /**< @brief Face marker for face consisting of second, third and fourth vertex */
    int nb;  /**< @brief Face marker for face consisting of first, third and fourth vertex */
    int nc;  /**< @brief Face marker for face consisting of first, second, and fourth vertex */
    int nd;  /**< @brief Face marker for face consisting of first, second, and third vertex */
};

class SurfaceMeshOld;

/** @brief Other data structure GemMesh (for volumetric meshes) */
class GemMesh {
public:

    // GemMesh constructors
    GemMesh() {}

    ~GemMesh() {}

    GemMesh(unsigned int num_vertices, unsigned int num_cells, bool higher_order = false);
    GemMesh(tetgenio& tetio);
    GemMesh(SurfaceMeshOld **surfmeshes, int num_meshes, char *tetgen_params);
    GemMesh(SurfaceMeshOld *surfmesh, char *tetgen_params);
    GemMesh(tetgenio *out, float radius, float centerx, float centery, float centerz, char *ActiveSite, int output_flag);

    // Methods working on a GemMesh
    void writeMcsf(char *filename);
    void writeDolfin(char *filename);
    void writeDiffpack(char *filename, int num_boundaries,
                       char *boundaries[]);
    void writeCarp(char *filename, int num_boundaries,
                   char *boundaries[]);
    void writeOFF(char *filename);

public:

    int  dim;              /**< @brief Geometrical dimension */
    int  dimii;            /**< @brief Topological dimension */
    int  num_vertices;     /**< @brief Number of vertices */
    int  num_cells;        /**< @brief Number of cells */
    bool higher_order;
    FETK_VX  *vv;          /**< @brief Vertex array */
    FETK_SS  *ss;          /**< @brief Cell array */
    INT6VECT *hs;          /**< @brief Higher order Cell array */
    SurfaceMeshOld *boundary; /**< @brief Boundary facets */
};


// Forward declarations
class tetgenio;

// GemMesh constructors
GemMesh* GemMesh_ctor(unsigned int num_vertices, unsigned int num_cells, bool higher_order = false);
GemMesh* GemMesh_fromTetgen(tetgenio& tetio);
GemMesh* GemMesh_fromSurfaceMeshOldes(SurfaceMeshOld **surfmeshes, int num_meshes, char *tetgen_params);
GemMesh* GemMesh_fromSurfaceMeshOld(SurfaceMeshOld *surfmesh, char *tetgen_params);
GemMesh* GemMesh_fromPdb(tetgenio *out, float radius, float centerx, float centery, float centerz,
                         char *ActiveSite, int output_flag);

// Methods working on a GemMesh
void GemMesh_writeMcsf(GemMesh *out, char *filename);
void GemMesh_writeDolfin(GemMesh *out, char *filename);
void GemMesh_writeDiffpack(GemMesh *out, char *filename, int num_boundaries, char *boundaries[]);
void GemMesh_writeCarp(GemMesh *out, char *filename, int num_boundaries, char *boundaries[]);
void GemMesh_dtor(GemMesh *);
void GemMesh_writeOFF(GemMesh *Gem_mesh, char *filename);


// void ReadActiveSiteFile(char*, unsigned int&, ATOM*&, unsigned int*&, float*&);
void  ReadActiveSiteFile(char *, unsigned int&, ATOM *&, unsigned int *&);
void  ReadRawiv(int *, int *, int *, float **, char *, float *, float *);
void  SurfaceExtract(TeTraMesh *, SurfaceMeshOld *);

#endif /* _BIOM_H_ */
