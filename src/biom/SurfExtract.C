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
 * File:     SurfExtract.C    < ... >
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Extract the surface mesh from a given tetrahedral mesh
 * ***************************************************************************
 */


#include <gamer/biom.h>
#include "gamercf.h"



/*
 * ***************************************************************************
 * Routine:  SurfaceExtract   < ... >
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Extract the surface mesh from a given tetrahedral mesh
 *           FIXME: It should be possible to speed up this...
 * ***************************************************************************
 */
void SurfaceExtract(TeTraMesh* volmesh, SurfaceMesh* surfmesh)
{
  int m,n,l;
  int start_index, face_num, num_vertices, vertex_index;
  int *vertex_indices;
  int a0, b0, c0, d0;
  int a, b, c, d;
  float x, y, z;
  char *boundary_flag;
  char surf_flag;
  unsigned int num_surfaces_vertices;
  INT3VECT* surfmesh_face; 

  // Initalize enough memory
  surfmesh_face = (INT3VECT*)malloc(sizeof(INT3VECT)*volmesh->num_cells);
 
  /* Classify the vertices into boundary or interior */
  boundary_flag = (char *)malloc(sizeof(char)*volmesh->num_vertices);
  for (n = 0; n < volmesh->num_vertices; n++) 
    boundary_flag[n] = 0;
  
  // Walk over all tetrahedrons of the Volume mesh
  face_num = 0;
  num_vertices = 0;
  for (n = 0; n < volmesh->num_cells; n++) {
    printf("%2.1f%% done \r", 100.0*(double)(n)/(double)(volmesh->num_cells-1));
    fflush(stdout);

    // Get the vertices of the tetrahedron
    a = volmesh->face[n].a;
    b = volmesh->face[n].b;
    c = volmesh->face[n].c;
    d = volmesh->face[n].d;
    
    // Walk over all tetrahedrons and check if the (a, b, c) face is shared by another tet
    surf_flag = 1;
    for (m = 0; m < volmesh->num_cells; m++) {
      if (m != n) {
	a0 = volmesh->face[m].a;
	b0 = volmesh->face[m].b;
	c0 = volmesh->face[m].c;
	d0 = volmesh->face[m].d;
	l = 0;
	if (a==a0 || a==b0 || a==c0 || a==d0)
	  l++;
	if (b==a0 || b==b0 || b==c0 || b==d0)
	  l++;
	if (c==a0 || c==b0 || c==c0 || c==d0)
	  l++;

	// m and n share a face
	if (l == 3) {
	  surf_flag = 0;
	  break;
	}
      }
    }
    if (surf_flag == 1) {
      num_vertices += (boundary_flag[a] == 0) ? 1 : 0;
      num_vertices += (boundary_flag[b] == 0) ? 1 : 0;
      num_vertices += (boundary_flag[c] == 0) ? 1 : 0;
      boundary_flag[a] = 1;
      boundary_flag[b] = 1;
      boundary_flag[c] = 1;
      surfmesh_face[face_num].a = a;
      surfmesh_face[face_num].b = c;
      surfmesh_face[face_num].c = b;
      face_num++;
    }

    // Walk over all tetrahedrons and check if the (a, b, d) face is shared by another tet
    surf_flag = 1;
    for (m = 0; m < volmesh->num_cells; m++) {
      if (m != n) {
	a0 = volmesh->face[m].a;
	b0 = volmesh->face[m].b;
	c0 = volmesh->face[m].c;
	d0 = volmesh->face[m].d;
	l = 0;
	if (a==a0 || a==b0 || a==c0 || a==d0)
	  l++;
	if (b==a0 || b==b0 || b==c0 || b==d0)
	  l++;
	if (d==a0 || d==b0 || d==c0 || d==d0)
	  l++;
	if (l == 3) {
	  surf_flag = 0;
	  break;
	}
      }
    }
    
    // m and n share a face
    if (surf_flag == 1) {
      num_vertices += (boundary_flag[a] == 0) ? 1 : 0;
      num_vertices += (boundary_flag[b] == 0) ? 1 : 0;
      num_vertices += (boundary_flag[d] == 0) ? 1 : 0;
      boundary_flag[a] = 1;
      boundary_flag[b] = 1;
      boundary_flag[d] = 1;
      surfmesh_face[face_num].a = a;
      surfmesh_face[face_num].b = b;
      surfmesh_face[face_num].c = d;
      face_num++;
    }

    // Walk over all tetrahedrons and check if the (a, c, d) face is shared by another tet
    surf_flag = 1;
    for (m = 0; m < volmesh->num_cells; m++) {
      if (m != n) {
	a0 = volmesh->face[m].a;
	b0 = volmesh->face[m].b;
	c0 = volmesh->face[m].c;
	d0 = volmesh->face[m].d;
	l = 0;
	if (a==a0 || a==b0 || a==c0 || a==d0)
	  l++;
	if (c==a0 || c==b0 || c==c0 || c==d0)
	  l++;
	if (d==a0 || d==b0 || d==c0 || d==d0)
	  l++;
	if (l == 3) {
	  surf_flag = 0;
	  break;
	}
      }
    }
    
    // m and n share a face
    if (surf_flag == 1) {
      num_vertices += (boundary_flag[a] == 0) ? 1 : 0;
      num_vertices += (boundary_flag[c] == 0) ? 1 : 0;
      num_vertices += (boundary_flag[d] == 0) ? 1 : 0;
      boundary_flag[a] = 1;
      boundary_flag[c] = 1;
      boundary_flag[d] = 1;
      surfmesh_face[face_num].a = a;
      surfmesh_face[face_num].b = d;
      surfmesh_face[face_num].c = c;
      face_num++;
    }

    // Walk over all tetrahedrons and check if the (b, c, d) face is shared by another tet
    surf_flag = 1;
    for (m = 0; m < volmesh->num_cells; m++) {
      if (m != n) {
	a0 = volmesh->face[m].a;
	b0 = volmesh->face[m].b;
	c0 = volmesh->face[m].c;
	d0 = volmesh->face[m].d;
	l = 0;
	if (b==a0 || b==b0 || b==c0 || b==d0)
	  l++;
	if (c==a0 || c==b0 || c==c0 || c==d0)
	  l++;
	if (d==a0 || d==b0 || d==c0 || d==d0)
	  l++;
	if (l == 3) {
	  surf_flag = 0;
	  break;
	}
      }
    }
    
    // m and n share a face
    if (surf_flag == 1) {
      num_vertices += (boundary_flag[b] == 0) ? 1 : 0;
      num_vertices += (boundary_flag[c] == 0) ? 1 : 0;
      num_vertices += (boundary_flag[d] == 0) ? 1 : 0;
      boundary_flag[b] = 1;
      boundary_flag[c] = 1;
      boundary_flag[d] = 1;
      surfmesh_face[face_num].a = b;
      surfmesh_face[face_num].b = c;
      surfmesh_face[face_num].c = d;
      face_num++;
    }
  }
  
  // Assign new memory to the surface mesh
  surfmesh->num_vertices = num_vertices;
  surfmesh->num_faces = face_num;
  surfmesh->vertex = (FLTVECT*)malloc(sizeof(FLTVECT)*num_vertices);
  surfmesh->face = (INT3VECT*)malloc(sizeof(INT3VECT)*face_num);
  
  // Assign vertex coordinates to surface mesh
  vertex_index = 0;
  vertex_indices = (int*)malloc(sizeof(int)*volmesh->num_vertices);
  for (n = 0; n < volmesh->num_vertices; n++) {
    if (boundary_flag[n]) {
      surfmesh->vertex[vertex_index].x = volmesh->vertex[n].x;
      surfmesh->vertex[vertex_index].y = volmesh->vertex[n].y;
      surfmesh->vertex[vertex_index].z = volmesh->vertex[n].z;
      vertex_indices[n] = vertex_index;
      vertex_index++;
    }
    else {
      vertex_indices[n] = -1;
    }
  }
  
  // Assign the correct face number to surface mesh
  for (n = 0; n < surfmesh->num_faces; n++) {
    surfmesh->face[n].a = vertex_indices[surfmesh_face[n].a];
    surfmesh->face[n].b = vertex_indices[surfmesh_face[n].b];
    surfmesh->face[n].c = vertex_indices[surfmesh_face[n].c];
  }
  
  printf("Surface Mesh Extracted: Nodes = %d, Faces = %d\n\n", surfmesh->num_vertices, surfmesh->num_faces);
      
  free(vertex_indices);
  free(boundary_flag);

  return;
}



/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_merge
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Merge two surface meshes into a single one
 * ***************************************************************************
 */
SurfaceMesh * SurfaceMesh_merge(SurfaceMesh *surfmesh_in0, SurfaceMesh *surfmesh_in1)
{
  SurfaceMesh *surfmesh;
  int i;
  
  surfmesh = SurfaceMesh_ctor(surfmesh_in0->num_vertices + surfmesh_in1->num_vertices, 
			      surfmesh_in0->num_faces + surfmesh_in1->num_faces);


  for (i = 0; i < surfmesh_in0->num_vertices; i++) {
    surfmesh->vertex[i].x = surfmesh_in0->vertex[i].x;
    surfmesh->vertex[i].y = surfmesh_in0->vertex[i].y;
    surfmesh->vertex[i].z = surfmesh_in0->vertex[i].z;
  }
  for (i = 0; i < surfmesh_in1->num_vertices; i++) {
    surfmesh->vertex[i+surfmesh_in0->num_vertices].x = surfmesh_in1->vertex[i].x;
    surfmesh->vertex[i+surfmesh_in0->num_vertices].y = surfmesh_in1->vertex[i].y;
    surfmesh->vertex[i+surfmesh_in0->num_vertices].z = surfmesh_in1->vertex[i].z;
  }
  
  for (i = 0; i < surfmesh_in0->num_faces; i++) {
    surfmesh->face[i].a = surfmesh_in0->face[i].a;
    surfmesh->face[i].b = surfmesh_in0->face[i].b;
    surfmesh->face[i].c = surfmesh_in0->face[i].c;
  }
  for (i = 0; i < surfmesh_in1->num_faces; i++) {
    surfmesh->face[i+surfmesh_in0->num_faces].a = surfmesh_in1->face[i].a+surfmesh_in0->num_vertices;
    surfmesh->face[i+surfmesh_in0->num_faces].b = surfmesh_in1->face[i].b+surfmesh_in0->num_vertices;
    surfmesh->face[i+surfmesh_in0->num_faces].c = surfmesh_in1->face[i].c+surfmesh_in0->num_vertices;
  }
  
  return surfmesh;

}

