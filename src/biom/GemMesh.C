/*
 * ***************************************************************************
 * GAMER = < Geometry-preserving Adaptive MeshER >
 * Copyright (C) 2007-2010 -- Michael Holst and Johan Hake
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
 * ***************************************************************************
 */

/* ***************************************************************************
 * File:     GemMesh.C    < ... >
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Create and destroy GemMesh data
 * ****************************************************************************
 */

#include <gamer/biom.h>
#include <gamer/tetgen.h>
#include "gamercf.h"
#include <vector>

// Declare internal GAMer methods
FLTVECT GetCrossProduct(SurfaceMesh *, int, int, int);

/*
 * ***************************************************************************
 * Routine:  same_face
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Check if the given numbers relate to the same face
 * ***************************************************************************
 */
bool same_face(int a, int b, int c, int* triface)
{
  int i;
  for (i=0; i<3; i++)
    if (a != triface[i] && b != triface[i] && c != triface[i])
      return false;
  return true;
}

/*
 * ***************************************************************************
 * Routine:  GemMesh_ctor
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Create memory for a GemMesh
 * ***************************************************************************
 */
GemMesh* GemMesh_ctor(unsigned int num_vertices, unsigned int num_cells, bool higher_order)
{
  GemMesh* gem_mesh = (GemMesh*)malloc(sizeof(GemMesh));
  gem_mesh->num_vertices = num_vertices;  
  gem_mesh->num_cells = num_cells; 
  gem_mesh->dim = 3;
  gem_mesh->dimii = 3;
  gem_mesh->vv = (FETK_VX *)malloc(sizeof(FETK_VX)*num_vertices);
  gem_mesh->ss = (FETK_SS *)malloc(sizeof(FETK_SS)*num_cells);
  gem_mesh->boundary = NULL;
  gem_mesh->higher_order = higher_order;
  
  if (higher_order)
    gem_mesh->hs = (INT6VECT *)malloc(sizeof(INT6VECT)*num_cells);
  else
    gem_mesh->hs = NULL;

  return gem_mesh;
}

/*
 * ***************************************************************************
 * Routine:  GemMesh_dtor
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Release memory from a GemMesh
 * ***************************************************************************
 */
void GemMesh_dtor(GemMesh* gem_mesh)
{
  free(gem_mesh->vv);
  free(gem_mesh->ss);
  free(gem_mesh);
  if (gem_mesh->boundary)
    SurfaceMesh_dtor(gem_mesh->boundary);
  
  if (gem_mesh->hs!=NULL)
    free(gem_mesh->hs);
  
}



/*
 * ***************************************************************************
 * Routine:  GemMesh_fromSurfaceMesh
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Use tetgen to generate a GemMesh from a surface mesh
 * ***************************************************************************
 */
GemMesh* GemMesh_fromSurfaceMesh(SurfaceMesh* surfmesh, char* tetgen_params)
{
  SurfaceMesh* surfmeshes[1] = {surfmesh};
  return GemMesh_fromSurfaceMeshes(surfmeshes, 1, tetgen_params);
}


/*
 * ***************************************************************************
 * Routine:  GemMesh_fromSurfaceMeshes
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Use tetgen to generate a GemMesh from a surface mesh
 * ***************************************************************************
 */
GemMesh* GemMesh_fromSurfaceMeshes(SurfaceMesh** surfmeshes, int num_surface_meshes, 
				   char* tetgen_params)
{
  
  unsigned int j, sm_index, num_vertices = 0, num_faces = 0;
  unsigned int num_regions = 0, num_holes = 0;
  SurfaceMesh* surfmesh;
  FLTVECT normal, midpoint;
  INT3VECT face0;

  tetgenio in, out;
  tetgenio::facet *f;
  tetgenio::polygon *p;
  in.firstnumber = 1;

  // Collect info from the surface meshes
  for (sm_index=0; sm_index<num_surface_meshes; sm_index++)
  {
    // Get the surface mesh
    surfmesh = surfmeshes[sm_index];

    // Check if neighborlist is created
    if (surfmesh->neighbor_list == NULL)
      SurfaceMesh_createNeighborlist(surfmesh);
    
    if (surfmesh->num_vertices == 0 || surfmesh->num_faces == 0)
    {
      printf("The %d:th surface mesh is empty. Abandoning "\
	     "tetrahedralization.\n", sm_index);
      return GemMesh_ctor(0, 0, false);
    }
      

    if (!surfmesh->closed)
    {
      printf("The %d:th surface mesh is not closed. Abandoning "\
	     "tetrahedralization.\n", sm_index);
      return GemMesh_ctor(0, 0, false);
    }
    
    // Bump counters
    if (!surfmesh->as_hole)
      num_regions += 1;
    else
      num_holes += 1;

    num_vertices += surfmesh->num_vertices;
    num_faces += surfmesh->num_faces;
  }

  // We need at least one region
  if (num_regions < 1)
  {
    printf("Expected at least one surface mesh which is not a hole.\n");
    return GemMesh_ctor(0, 0, false);
  }
  
  // Assign memory
  printf("num vertices: %d\n", num_vertices);
  printf("num faces: %d\n", num_faces);
  in.numberofpoints = num_vertices;
  in.pointlist = new REAL[in.numberofpoints * 3];
  
  in.numberoffacets = num_faces;
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];
  
  printf("num regions: %d\n", num_regions);
  in.numberofregions = num_regions;
  in.regionlist = new REAL[num_regions*5];

  if (num_holes > 0)
  {
    printf("num holes: %d\n", num_holes);
    in.numberofholes = num_holes;
    in.holelist = new REAL[num_holes*3];
  }

  // Reset counters
  num_vertices = num_faces = num_regions = num_holes = 0;
  
  // Iterate over the surface meshes and assign info
  for (sm_index=0; sm_index<num_surface_meshes; sm_index++)
  {
    
    // Get the surface mesh
    surfmesh = surfmeshes[sm_index];
    
    // Assign vertex information
    for (j = 0; j < surfmesh->num_vertices; j++) 
    {
      in.pointlist[j*3 + 0 + num_vertices*3] = surfmesh->vertex[j].x;
      in.pointlist[j*3 + 1 + num_vertices*3] = surfmesh->vertex[j].y;
      in.pointlist[j*3 + 2 + num_vertices*3] = surfmesh->vertex[j].z;
    }
  
    // Assign face information
    for (j = 0; j < surfmesh->num_faces; j++) 
    {
      f = &in.facetlist[j+num_faces];
      f->holelist = (REAL *) NULL;
      f->numberofholes = 0;
      f->numberofpolygons = 1;
      f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
      p = &f->polygonlist[0];
      p->numberofvertices = 3;
      p->vertexlist = new int[p->numberofvertices];
      p->vertexlist[0] = surfmesh->face[j].a + in.firstnumber + num_vertices;
      p->vertexlist[1] = surfmesh->face[j].b + in.firstnumber + num_vertices;
      p->vertexlist[2] = surfmesh->face[j].c + in.firstnumber + num_vertices;
      in.facetmarkerlist[j+num_faces] = surfmesh->face[j].m == -1 ? 0 : surfmesh->face[j].m;
    }

    // Get region or hole coordinate by using a face
    face0 = surfmesh->face[0];

    // Face normal 
    normal = GetCrossProduct(surfmesh, face0.a, face0.b, face0.c);
    
    // Midpoint of face
    midpoint.x = surfmesh->vertex[face0.a].x/3;
    midpoint.y = surfmesh->vertex[face0.a].y/3;
    midpoint.z = surfmesh->vertex[face0.a].z/3;

    midpoint.x += surfmesh->vertex[face0.b].x/3;
    midpoint.y += surfmesh->vertex[face0.b].y/3;
    midpoint.z += surfmesh->vertex[face0.b].z/3;

    midpoint.x += surfmesh->vertex[face0.c].x/3;
    midpoint.y += surfmesh->vertex[face0.c].y/3;
    midpoint.z += surfmesh->vertex[face0.c].z/3;

    // Weight based on length
    double weight = sqrt(pow(surfmesh->vertex[face0.a].x-\
			     surfmesh->vertex[face0.b].x, 2) + \
			 pow(surfmesh->vertex[face0.a].y-\
			     surfmesh->vertex[face0.b].y, 2) + \
			 pow(surfmesh->vertex[face0.a].z-\
			     surfmesh->vertex[face0.b].z, 2));

    // Add region information
    if (!surfmesh->as_hole)
    {
      // Calculate region info using face information
      // FIXME: Why should we add the normal values?
      // FIXME: If normal is pointing outwards we should substract to 
      // FIXME: get a point inside the surface mesh....
      in.regionlist[num_regions*5 + 0] = midpoint.x + normal.x*weight;
      in.regionlist[num_regions*5 + 1] = midpoint.y + normal.y*weight;
      in.regionlist[num_regions*5 + 2] = midpoint.z + normal.z*weight;
      
      // Set region marker (attribute)
      in.regionlist[num_regions*5 + 3] = surfmesh->marker;

      // If volume constraints
      if (surfmesh->use_volume_constraint)
      {
	printf("Volume constraint: %.5f\n", surfmesh->volume_constraint);
	in.regionlist[num_regions*5 + 4] = surfmesh->volume_constraint;
      }
      else
	in.regionlist[num_regions*5 + 4] = -1;
      
      // Bump counter
      num_regions += 1;
    }
    else
    {
      // Calculate hole info using face information
      in.holelist[num_holes*3 + 0] = midpoint.x + normal.x*weight;
      in.holelist[num_holes*3 + 1] = midpoint.y + normal.y*weight;
      in.holelist[num_holes*3 + 2] = midpoint.z + normal.z*weight;

      // Bump counter
      num_holes += 1;
    }
    
    // Bump face and vertex counters
    num_vertices += surfmesh->num_vertices;
    num_faces += surfmesh->num_faces;

  }  
  
  // Add boundary marker on each node
  // FIXME: Why? Aren't the markers set on the generatedmesh?
  in.pointmarkerlist = new int[in.numberofpoints];
  for (j = 0; j < in.numberofpoints; j++)
    in.pointmarkerlist[j] = 1;

  // Debug
  in.save_nodes("plc");
  in.save_poly("plc");

  // Call TetGen
  tetrahedralize(tetgen_params, &in, &out, NULL);

  // Debug
  out.save_nodes("result");
  out.save_elements("result");
  out.save_faces("result");

  // Convert to a generalized tetrahedral mesh structure
  return GemMesh_fromTetgen(out);

}


/*
 * ***************************************************************************
 * Routine:  GemMesh_fromTetgen
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Generate a GemMesh structure 
 * ***************************************************************************
 */
/*
FIXME: Are we using the neighbor list at all?
 */
GemMesh* GemMesh_fromTetgen(tetgenio& tetio)
{
  GemMesh *gem_mesh;
  unsigned int i, j, k;
  FILE *fout;
  long node_index;
  int *face_type, tetra_node[4];
  float x,y,z;
  int cell_type, boundary_face=0;
  float dist;
  long a,b,c,d;
  float ax,ay,az;
  float bx,by,bz;
  float cx,cy,cz;
  float dx,dy,dz;
  int active, exterior_cell_type;
  long neighbor, neightype;
  long *map_w_o,*map_w_i,*map_w_t; 
  Triface** tet_to_triface;
  int tet0, tet1, tet2;
  int *tetp;
  bool on_boundary[2];
  int* vertex_markers;
  Triface* triface;
  
  // Check for higher order cells
  const bool higher_order = tetio.numberofcorners==10;

  // Initialization and memory allocation
  gem_mesh = GemMesh_ctor(tetio.numberofpoints, tetio.numberoftetrahedra, higher_order);
  tet_to_triface = new Triface*[tetio.numberoftetrahedra];
  vertex_markers = new int[tetio.numberofpoints];

  for (i = 0; i < tetio.numberoftetrahedra; i++)
    tet_to_triface[i] = NULL;
  
  for (i = 0; i < tetio.numberofpoints; i++)
    vertex_markers[i] = 0;
  
  // Using the markers from the trifaces to mark vertices beeing on the boundary
  SurfaceMesh* boundary = SurfaceMesh_ctor(0, tetio.numberoftrifaces);
  for (i = 0; i < tetio.numberoftrifaces; i++)
  {
    tetp = &tetio.adjtetlist[i*2];
    on_boundary[0] = tetp[0] >= 0 && tetp[1] < 0;
    on_boundary[1] = tetp[1] >= 0 && tetp[0] < 0;
    
    // Iterate over the adjacent tets and register stored face information
    for (j=0; j<2; j++)
    {
      // If no tet continue
      if (tetp[j] < 0)
	continue;
      
      triface = new Triface;
      triface->next = NULL;
      triface->points[0] = tetio.trifacelist[i*3]-1;
      triface->points[1] = tetio.trifacelist[i*3+1]-1;
      triface->points[2] = tetio.trifacelist[i*3+2]-1;
      
      // If we are on a boundary face register it
      if (on_boundary[j])
      {
	boundary->face[boundary_face].a = triface->points[0];
	boundary->face[boundary_face].b = triface->points[1];
	boundary->face[boundary_face].c = triface->points[2];
	if (tetio.trifacemarkerlist != NULL)
	  boundary->face[boundary_face].m = tetio.trifacemarkerlist[i];
	else
	  boundary->face[boundary_face].m = 1;
	boundary_face++;
      }
      
      if (tetio.trifacemarkerlist != NULL)
	triface->marker = tetio.trifacemarkerlist[i];
      else
	triface->marker = 1;
      
      triface->next = tet_to_triface[tetp[j]-1];
      tet_to_triface[tetp[j]-1] = triface;
      
      vertex_markers[triface->points[0]] = triface->marker;
      vertex_markers[triface->points[1]] = triface->marker;
      vertex_markers[triface->points[2]] = triface->marker;
    }
  }
  
  // Store the actuall number of boundary faces
  boundary->num_faces = boundary_face;

  // Generate facetype information
  face_type = new int[tetio.numberoftetrahedra*4];

  for (i = 0; i < tetio.numberoftetrahedra; i++) {
    
    // Zero out all face types
    for (j = 0; j < 4; j++) 
      face_type[i*4+j] = 0;

    // Check tetra orientation
    a = tetio.tetrahedronlist[i * tetio.numberofcorners + 0] - 1;
    b = tetio.tetrahedronlist[i * tetio.numberofcorners + 1] - 1;
    c = tetio.tetrahedronlist[i * tetio.numberofcorners + 2] - 1;
    d = tetio.tetrahedronlist[i * tetio.numberofcorners + 3] - 1;
    ax = tetio.pointlist[a * 3];
    ay = tetio.pointlist[a * 3 + 1];
    az = tetio.pointlist[a * 3 + 2];
    bx = tetio.pointlist[b * 3];
    by = tetio.pointlist[b * 3 + 1];
    bz = tetio.pointlist[b * 3 + 2];
    cx = tetio.pointlist[c * 3];
    cy = tetio.pointlist[c * 3 + 1];
    cz = tetio.pointlist[c * 3 + 2];
    dx = tetio.pointlist[d * 3];
    dy = tetio.pointlist[d * 3 + 1];
    dz = tetio.pointlist[d * 3 + 2];
    bx -= ax;
    by -= ay;
    bz -= az;
    cx -= ax;
    cy -= ay;
    cz -= az;
    dx -= ax;
    dy -= ay;
    dz -= az;
    ax = by*cz-bz*cy;
    ay = bz*cx-bx*cz;
    az = bx*cy-by*cx;
    
    if (ax*dx+ay*dy+az*dz < 0 && !higher_order)
    {
      tetio.tetrahedronlist[i * 4 + 1] = c;
      tetio.tetrahedronlist[i * 4 + 2] = b;
      k = tetio.neighborlist[i * 4 + 1];
      tetio.neighborlist[i * 4 + 1] = tetio.neighborlist[i * 4 + 2];
      tetio.neighborlist[i * 4 + 2] = k;
    }

    // Check if we have a marked face
    while (tet_to_triface[i] != NULL)
    {

      // If so find out which face is on the boundary
      if (same_face(b, c, d, tet_to_triface[i]->points))
      {
	face_type[i*4] = tet_to_triface[i]->marker;
	triface = tet_to_triface[i];
	tet_to_triface[i] = triface->next;
	delete triface;
	continue;
      }

      if (same_face(a, c, d, tet_to_triface[i]->points))
      {
	face_type[i*4+1] = tet_to_triface[i]->marker;
	triface = tet_to_triface[i];
	tet_to_triface[i] = triface->next;
	delete triface;
	continue;
      }

      if (same_face(a, b, d, tet_to_triface[i]->points))
      {
	face_type[i*4+2] = tet_to_triface[i]->marker;
	triface = tet_to_triface[i];
	tet_to_triface[i] = triface->next;
	delete triface;
	continue;
      }

      if (same_face(a, b, c, tet_to_triface[i]->points))
      {
	face_type[i*4+3] = tet_to_triface[i]->marker;
	triface = tet_to_triface[i];
	tet_to_triface[i] = triface->next;
	delete triface;
	continue;
      }
    }
  }

  // Start output mesh
  gem_mesh->num_vertices = tetio.numberofpoints;
  gem_mesh->num_cells = tetio.numberoftetrahedra;

  // Write all nodes to a GemMesh structure
  for (i = 0; i < tetio.numberofpoints; i++)
  {
      gem_mesh->vv[i].id = i;
      gem_mesh->vv[i].chrt = vertex_markers[i]; // This might be ambigious
      gem_mesh->vv[i].x = tetio.pointlist[i * 3];
      gem_mesh->vv[i].y = tetio.pointlist[i * 3 + 1];
      gem_mesh->vv[i].z = tetio.pointlist[i * 3 + 2];
  }

  // Write all tets to a GemMesh structure
  for (i = 0; i < tetio.numberoftetrahedra; i++) 
  {
    gem_mesh->ss[i].id = i;
    gem_mesh->ss[i].grp = 0;
  
    if (tetio.numberoftetrahedronattributes > 0)
      gem_mesh->ss[i].mat = (int)tetio.tetrahedronattributelist[i*tetio.numberoftetrahedronattributes];
    else
      gem_mesh->ss[i].mat = 0;

    gem_mesh->ss[i].fa = face_type[4*i+0];
    gem_mesh->ss[i].fb = face_type[4*i+1];
    gem_mesh->ss[i].fc = face_type[4*i+2];
    gem_mesh->ss[i].fd = face_type[4*i+3];

    // NOTE: Index in tetrahedronlist starts from 1
    gem_mesh->ss[i].na = tetio.tetrahedronlist[i*tetio.numberofcorners+0] - 1;
    gem_mesh->ss[i].nb = tetio.tetrahedronlist[i*tetio.numberofcorners+1] - 1;
    gem_mesh->ss[i].nc = tetio.tetrahedronlist[i*tetio.numberofcorners+2] - 1;
    gem_mesh->ss[i].nd = tetio.tetrahedronlist[i*tetio.numberofcorners+3] - 1;

    // Add higher order cell info
    if (higher_order)
    {
      gem_mesh->hs[i].a = tetio.tetrahedronlist[i*tetio.numberofcorners+4] - 1;
      gem_mesh->hs[i].b = tetio.tetrahedronlist[i*tetio.numberofcorners+5] - 1;
      gem_mesh->hs[i].c = tetio.tetrahedronlist[i*tetio.numberofcorners+6] - 1;
      gem_mesh->hs[i].d = tetio.tetrahedronlist[i*tetio.numberofcorners+7] - 1;
      gem_mesh->hs[i].e = tetio.tetrahedronlist[i*tetio.numberofcorners+8] - 1;
      gem_mesh->hs[i].f = tetio.tetrahedronlist[i*tetio.numberofcorners+9] - 1;
    }
    
  }

  // Store boundary
  gem_mesh->boundary = boundary;

  // Release memory
  delete[] face_type;
  delete[] tet_to_triface;
  delete[] vertex_markers;

  return gem_mesh;
}

