/*
 * ***************************************************************************
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
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     SurfRefine.C
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Refine the surface mesh by subdividing each triangle into four
 * ***************************************************************************
 */
            
#include <gamer/biom.h>
#include "gamercf.h"


/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_refine
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Refine a surface mesh and recalculate the neighbors
 * ***************************************************************************
 */
void SurfaceMesh_refine(SurfaceMesh* surfmesh)
{
  unsigned int local_num_edges = 0, total_num_edges = 0, edge_num = 0;
  unsigned int min_vertex_num, max_vertex_num, face_num;
  unsigned int k, m, n, nv;
  unsigned int a, b, c;
  unsigned int* num_edges;
  unsigned int* offsets;
  unsigned int* vertex2edge;
  unsigned int local_vertices[3], local_additional_vertices[3];
  NPNT3* ngr;
  SurfaceMesh* surfmesh_refine;
  float ax, ay, az;
  float nx, ny, nz;
  
  // Check if neighborlist is created
  if (surfmesh->neighbor_list == NULL)
    SurfaceMesh_createNeighborlist(surfmesh);

  NPNT3** neighbor_list = surfmesh->neighbor_list;

  // Assign the number of vertices in the original mesh
  nv = surfmesh->num_vertices;

  // Create an array with the number of edges associated with each vertex
  num_edges = (unsigned int*)malloc(sizeof(unsigned int)*nv);

  // Create an array with the offsets into the vertex2edge array for each vertex
  offsets = (unsigned int*)malloc(sizeof(unsigned int)*nv);

  // Iterate over all vertices and collect edges
  for (n = 0; n < nv; n++){
    offsets[n] = total_num_edges;
    local_num_edges = 0;
    ngr = neighbor_list[n];
    while (ngr != NULL) {
      // If n is smaller than ngr->a we have an edge
      if (n < ngr->a){
	total_num_edges++;
	local_num_edges++;
      }
      ngr = ngr->next;
    }
    num_edges[n] = local_num_edges;
    
  }
  
  // Create memory for the refined mesh
  surfmesh_refine = SurfaceMesh_ctor(nv + total_num_edges, surfmesh->num_faces*4);
  surfmesh_refine->num_vertices = nv;
  surfmesh_refine->num_faces = surfmesh->num_faces;
  
  // Copy the original mesh to the new mesh
  for (n = 0; n < nv; n++) {
    surfmesh_refine->vertex[n].x = surfmesh->vertex[n].x;
    surfmesh_refine->vertex[n].y = surfmesh->vertex[n].y;
    surfmesh_refine->vertex[n].z = surfmesh->vertex[n].z;
  }
  
  for (n = 0; n < surfmesh->num_faces; n++) {
    surfmesh_refine->face[n].a = surfmesh->face[n].a;
    surfmesh_refine->face[n].b = surfmesh->face[n].b;
    surfmesh_refine->face[n].c = surfmesh->face[n].c;
  }
  
  // Create the map from vertices to edges
  vertex2edge = (unsigned int*)malloc(sizeof(unsigned int)*total_num_edges);

  // Iterate over all vertices and split edges
  for (n = 0; n < nv; n++){

    // Get the coordinates of vertex n
    nx = surfmesh_refine->vertex[n].x;
    ny = surfmesh_refine->vertex[n].y;
    nz = surfmesh_refine->vertex[n].z;

    ngr = neighbor_list[n];
    while (ngr != NULL) {
   
      // If n is smaller than ngr->a we have an edge
      if (n < ngr->a){

	// Add the value of the opposite vertex to the map
	vertex2edge[edge_num] = ngr->a;
	
	// Get the coordinates of vertex ngr->a
	ax = surfmesh_refine->vertex[ngr->a].x;
	ay = surfmesh_refine->vertex[ngr->a].y;
	az = surfmesh_refine->vertex[ngr->a].z;
	
	// Add the new vertex coordinates of the splitted edge
	surfmesh_refine->vertex[nv + edge_num].x = 0.5*(ax + nx);
	surfmesh_refine->vertex[nv + edge_num].y = 0.5*(ay + ny);
	surfmesh_refine->vertex[nv + edge_num].z = 0.5*(az + nz);
	
	// Increase the edge number
	edge_num++;
      }
      ngr = ngr->next;
    }
  }

  // A counter for adding new faces
  face_num = surfmesh_refine->num_faces;

  // Iterate over faces and add information of the refined face
  for (n = 0; n < surfmesh_refine->num_faces; n++) {
    local_vertices[0] = surfmesh_refine->face[n].a;
    local_vertices[1] = surfmesh_refine->face[n].b;
    local_vertices[2] = surfmesh_refine->face[n].c;

    // Iterate over the vertices and find the edges
    for (m = 0; m < 3; m++){
      min_vertex_num = min(local_vertices[m], local_vertices[(m+1)%3]);
      max_vertex_num = max(local_vertices[m], local_vertices[(m+1)%3]);
      
      // Find the edge number that fit the pair of vertices
      for (k = 0; k < num_edges[min_vertex_num]; k++)
	if (vertex2edge[offsets[min_vertex_num] + k] == max_vertex_num)
	  break;
      
      // The edge number represents the number of the added vertex plus the 
      // number of original vertices
      local_additional_vertices[m] = nv + offsets[min_vertex_num] + k;

    }
    
    // Add information of the four new faces

    // First the mid face
    surfmesh_refine->face[n].a = local_additional_vertices[0];
    surfmesh_refine->face[n].b = local_additional_vertices[1];
    surfmesh_refine->face[n].c = local_additional_vertices[2];

    // Then the three corner faces
    for (m = 0; m < 3; m++){
      surfmesh_refine->face[face_num].a = local_vertices[m];
      surfmesh_refine->face[face_num].b = local_additional_vertices[m];
      surfmesh_refine->face[face_num].c = local_additional_vertices[(m+2)%3];
      face_num++;
    }
  }

  // Release memory
  free(num_edges);
  free(offsets);
  free(vertex2edge);

  // Update number information
  surfmesh_refine->num_vertices += total_num_edges;
  surfmesh_refine->num_faces *= 4;
  
  // Release old data
  SurfaceMesh_releaseData(surfmesh);

  // Assign the refined mesh to the passed
  surfmesh->num_vertices = surfmesh_refine->num_vertices;
  surfmesh->num_faces = surfmesh_refine->num_faces;
  surfmesh->vertex = surfmesh_refine->vertex;
  surfmesh->face = surfmesh_refine->face;
  
  // Free memory of refined surface mesh struct
  free(surfmesh_refine);

  // Recreate the neigborlist
  SurfaceMesh_createNeighborlist(surfmesh);

}
