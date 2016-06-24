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
 * File:     SurfaceMesh.C    < ... >
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Create and destroy SurfaceMesh data
 * ****************************************************************************
 */

#include <gamer/biom.h>
#include <vector>
#include <algorithm>
#include <utility>
#include "gamercf.h"
#include <cmath>

// Declare internal GAMer methods
EIGENVECT GetEigenVector(SurfaceMesh *, int, FLTVECT *,float *);

/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_ctor
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Create a surface mesh instance
 * ***************************************************************************
 */
SurfaceMesh* SurfaceMesh_ctor(unsigned int num_vertices, unsigned int num_faces)
{
  // Allocate memory and initialize variables
  SurfaceMesh* surfmesh = (SurfaceMesh*)malloc(sizeof(SurfaceMesh));
  surfmesh->num_vertices = num_vertices;
  surfmesh->num_faces = num_faces;
  if (num_vertices)
    surfmesh->vertex = (FLTVECT*)malloc(sizeof(FLTVECT)*num_vertices);
  else
    surfmesh->vertex = NULL;
  if (num_faces)
    surfmesh->face = (INT3VECT*)malloc(sizeof(INT3VECT)*num_faces);
  else
    surfmesh->face = NULL;
  surfmesh->neighbor = NULL;
  surfmesh->neighbor_list = NULL;
  surfmesh->avglen = 0.;
  surfmesh->min[0] = 0.; surfmesh->min[1] = 0.; surfmesh->min[2] = 0.;
  surfmesh->max[0] = 0.; surfmesh->max[1] = 0.; surfmesh->max[2] = 0.;

  // Initialize SurfaceMesh structures
  for (int n = 0; n < num_vertices; n++)
  {
    FLTVECT& vert = surfmesh->vertex[n];
    vert.x = vert.y = vert.z = vert.m = 0;
    vert.sel = true;
  }
  
  for (int n = 0; n < num_faces; n++)
  {
    INT3VECT& face = surfmesh->face[n];
    face.a = face.b = face.c = face.m = 0;
    face.sel = true;
  }

  // Initialize domain data
  surfmesh->closed = true;
  surfmesh->marker = 1;
  surfmesh->volume_constraint = 100; 
  surfmesh->use_volume_constraint = false; 
  surfmesh->as_hole = false; 

  return surfmesh;
}

/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_dtor
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Release memory from a SurfaceMesh instance
 * ***************************************************************************
 */
void SurfaceMesh_dtor(SurfaceMesh* surfmesh)
{
  // Relase data memory
  SurfaceMesh_releaseData(surfmesh);

  // Release memory for the SurfaceMesh struct
  free(surfmesh);
}

/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_releaseData
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Release data memory from a SurfaceMesh instance
 * ***************************************************************************
 */
void SurfaceMesh_releaseData(SurfaceMesh* surfmesh)
{
  // Free allocated memory
  if (surfmesh->vertex)
    free(surfmesh->vertex);
  if (surfmesh->face)
    free(surfmesh->face);

  // Destroy neighbor_list
  SurfaceMesh_destroyNeighborlist(surfmesh);
  
}

/*
 * ***************************************************************************
 * SubRoutine:  SurfaceMesh_createNeighborList
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com) and Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Create a neighbor list from a surface mesh
 * ***************************************************************************
 */
void SurfaceMesh_createNeighborlist(SurfaceMesh* surfmesh)
{
  
  int m, n, a0, b0;
  int a, b, c, d;
  NPNT3 *first_ngr, *second_ngr, *tmp_ngr, *last_ngr;
  NPNT3 **neighbor_list;
  bool closed;

  // Destroy any exsisting neighborlist
  SurfaceMesh_destroyNeighborlist(surfmesh);

  // Create an array of NPNT3, used to store 
  neighbor_list = (NPNT3 **)malloc(sizeof(NPNT3 *)*surfmesh->num_vertices);

  // Initialize the neighbor list
  for (n = 0; n < surfmesh->num_vertices; n++)
  {
    neighbor_list[n] = NULL;

    // Default mark all vertices for deletion
    surfmesh->vertex[n].m = -1;
  }

  // Iterate over the faces and collect line segments (a, b) and its connection 
  // to a face (c). Save the line segment so it forms a counter clockwise triangle 
  // with the origin vertex 
  int num_connected = 0;
  for (n = 0; n < surfmesh->num_faces; n++) 
  {
  
    a = surfmesh->face[n].a;
    b = surfmesh->face[n].b;
    c = surfmesh->face[n].c;
   
    if (a == b || b == c || a == c)
      printf("Face %d include vertices with same indices (%d, %d, %d).\n", n, a, b, c);

    first_ngr = (NPNT3 *)malloc(sizeof(NPNT3));
    first_ngr->a = b;
    first_ngr->b = c;
    first_ngr->c = n;
    first_ngr->next = neighbor_list[a];
    neighbor_list[a] = first_ngr;

    // Mark vertex as connected
    if (surfmesh->vertex[a].m<0)
    {
      surfmesh->vertex[a].m = 0;
      num_connected += 1;
    }

    first_ngr = (NPNT3 *)malloc(sizeof(NPNT3));
    first_ngr->a = c;
    first_ngr->b = a;
    first_ngr->c = n;
    first_ngr->next = neighbor_list[b];
    neighbor_list[b] = first_ngr;

    // Mark vertex as connected
    if (surfmesh->vertex[b].m<0)
    {
      surfmesh->vertex[b].m = 0;
      num_connected += 1;
    }
    
    first_ngr = (NPNT3 *)malloc(sizeof(NPNT3));
    first_ngr->a = a;
    first_ngr->b = b;
    first_ngr->c = n;
    first_ngr->next = neighbor_list[c];
    neighbor_list[c] = first_ngr;

    // Mark vertex as connected
    if (surfmesh->vertex[c].m<0)
    {
      surfmesh->vertex[c].m = 0;
      num_connected += 1;
    }
    
  }

  // Check if there are vertices which are not connect to any face
  if (num_connected<surfmesh->num_vertices)
  {
    // Attach the neighborlist to the surfmesh and destroy it
    surfmesh->neighbor_list = neighbor_list;
    SurfaceMesh_destroyNeighborlist(surfmesh);
    
    // Remove unconnected vertices
    SurfaceMesh_removeUnconnectedVertices(surfmesh);
    
    // Re-create neighbors
    SurfaceMesh_createNeighborlist(surfmesh);
    
    return;
    
  }

  // Order the neighbors so they are connected counter clockwise
  for (n = 0; n < surfmesh->num_vertices; n++) 
  {
    a0 = b0 = -1;
    first_ngr = neighbor_list[n];
    c = first_ngr->a;
    d = first_ngr->b;
    while (first_ngr != NULL) 
    {
      a = first_ngr->a;
      b = first_ngr->b;
      
      second_ngr = first_ngr->next;
      while (second_ngr != NULL) 
      {
	a0 = second_ngr->a;
	b0 = second_ngr->b;
	if (a0==b && b0!=a) 
	{
	  tmp_ngr = first_ngr;
	  while (tmp_ngr != NULL) 
	  {
	    if (tmp_ngr->next == second_ngr) 
	    {
	      tmp_ngr->next = second_ngr->next;
	      break;
	    }
	    tmp_ngr = tmp_ngr->next;
	  }
	  tmp_ngr = first_ngr->next;
 	  first_ngr->next = second_ngr;
	  second_ngr->next = tmp_ngr;
	  break;
	}

	second_ngr = second_ngr->next;

      }
      
      first_ngr = first_ngr->next;
      
    }
    
    // Check that the neighbor list is connected
    tmp_ngr = neighbor_list[n];

    closed = true;
    while (tmp_ngr->next != NULL)
    {
      // Check that we are connected
      if (tmp_ngr->b != tmp_ngr->next->a)
      {
	if (closed)
	{
	  printf("Polygons connected to vertex %d are not closed (interupted):"
		 " (%.2f, %.2f, %.2f)\n", n, surfmesh->vertex[n].x,
		 surfmesh->vertex[n].y, surfmesh->vertex[n].z);
	}
	  
	// Do not bail, just register the vertex to not be done anything with
	surfmesh->vertex[n].sel = false;

	closed = false;
      }
      
      // Step one face forward
      tmp_ngr = tmp_ngr->next;
    }

    // Check if the list forms a closed ring
    if (closed && b0 != c)
    {
      printf("Polygons connected to vertex %d are not closed (not closed):"
	     " (%.2f, %.2f, %.2f)\n", n, surfmesh->vertex[n].x,
	     surfmesh->vertex[n].y, surfmesh->vertex[n].z);
      
      // Do not bail, just register the vertex to not be done anything with
      surfmesh->vertex[n].sel = false;
      
      closed = false;
    }
    
    if (!closed)
    {
      surfmesh->closed = false;

    }
    
  }

  // Attach the neighborlist to the surfmesh
  surfmesh->neighbor_list = neighbor_list;
}

/*
 * ***************************************************************************
 * SubRoutine:  SurfaceMesh_destroyNeighborlist
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Release memory from a neighborlist
 * ***************************************************************************
 */
void SurfaceMesh_destroyNeighborlist(SurfaceMesh* surfmesh)
{
  unsigned int n;
  NPNT3* first_ngr = NULL;
  NPNT3* tmp_ngr = NULL;

  if (surfmesh->neighbor_list != NULL)
  {
    // Release the single neighbors 
    for (n = 0; n < surfmesh->num_vertices; n++) 
    {
      first_ngr = surfmesh->neighbor_list[n];
      while (first_ngr != NULL) 
      {
	tmp_ngr = first_ngr->next;
	free(first_ngr);
	first_ngr = tmp_ngr;
      }
    } 

    // Free the array of pointers
    free(surfmesh->neighbor_list);
    surfmesh->neighbor_list = NULL;
  }
}

/*
 * ***************************************************************************
 * SubRoutine:  SurfaceMesh_flipNormals
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Release memory from a neighborlist
 * ***************************************************************************
 */
void SurfaceMesh_flipNormals(SurfaceMesh* surfmesh)
{

  // First destroy neighborlist
  SurfaceMesh_destroyNeighborlist(surfmesh);

  // Iterate over faces and flip the vertices
  for (int n = 0; n < surfmesh->num_faces; n++) 
  {
    const int tmp = surfmesh->face[n].a;
    surfmesh->face[n].a = surfmesh->face[n].b;
    surfmesh->face[n].b = tmp;
  }

}


/*
 * ***************************************************************************
 * SubRoutine: SurfaceMesh_remove_unconnected_vertices
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Remove vertices which are not connected to any face
 * ***************************************************************************
 */
void SurfaceMesh_removeUnconnectedVertices(SurfaceMesh* surfmesh)
{
  // Collect statistics
  int num_to_be_removed = 0 ;
  std::vector<int> to_be_removed(surfmesh->num_vertices);
  for (int n = 0; n < surfmesh->num_vertices; n++)
  {
    if (surfmesh->vertex[n].m<0)
      num_to_be_removed++;
    to_be_removed[n] = num_to_be_removed;
  }
  
  printf("Removing %d vertices.\n", num_to_be_removed);

  // Move vertices forward
  for (int n = 0; n < surfmesh->num_vertices; n++)
  {
    // If a vertex is to be removed
    if ((n==0 && to_be_removed[n] != 0) || 
	(n!=0 && to_be_removed[n-1]!=to_be_removed[n]) )
      continue;
    
    // Move vertices forward
    surfmesh->vertex[n-to_be_removed[n]].x = surfmesh->vertex[n].x;
    surfmesh->vertex[n-to_be_removed[n]].y = surfmesh->vertex[n].y;
    surfmesh->vertex[n-to_be_removed[n]].z = surfmesh->vertex[n].z;
    surfmesh->vertex[n-to_be_removed[n]].sel = surfmesh->vertex[n].sel;
    surfmesh->vertex[n-to_be_removed[n]].m = surfmesh->vertex[n].m;
  }
    
  // Fix face offset
  for (int n = 0; n < surfmesh->num_faces; n++)
  {
    surfmesh->face[n].a = surfmesh->face[n].a - to_be_removed[surfmesh->face[n].a];
    surfmesh->face[n].b = surfmesh->face[n].b - to_be_removed[surfmesh->face[n].b];
    surfmesh->face[n].c = surfmesh->face[n].c - to_be_removed[surfmesh->face[n].c];
  }
    
  // Adjust num_vertices
  surfmesh->num_vertices -= num_to_be_removed;

}

/*
 * ***************************************************************************
 * SubRoutine: SurfaceMesh_deleteVertices
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Remove vertices which are marked with -1 as a vertex marker
 * ***************************************************************************
 */
void SurfaceMesh_deleteVertices(SurfaceMesh* surfmesh)
{
  // Mark faces connected to vertices for deletion
  for (int n=0; n < surfmesh->num_faces; n++)
  {
    if (surfmesh->vertex[surfmesh->face[n].a].m < 0 ||
	surfmesh->vertex[surfmesh->face[n].b].m < 0 ||
	surfmesh->vertex[surfmesh->face[n].c].m < 0 )
    {
      surfmesh->face[n].m = -1;
    }
  }

  // Delete marked faces
  SurfaceMesh_deleteFaces(surfmesh);

}

/*
 * ***************************************************************************
 * SubRoutine: SurfaceMesh_deleteFaces
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Remove faces which are marked with -1 as a face marker
 * ***************************************************************************
 */
void SurfaceMesh_deleteFaces(SurfaceMesh* surfmesh)
{

  // Iterate over vertices and mark all for deletion
  for (int n=0; n < surfmesh->num_vertices; n++)
    surfmesh->vertex[n].m = -1;
  
  // Delete faces connected to vertices
  int num_removed = 0;
  for (int n=0; n < surfmesh->num_faces; n++)
  {
    // Check for removal of face
    if (surfmesh->face[n].m < 0)
      num_removed += 1;
    else 
    {
      // If any previous face has been marked for deletion
      if (num_removed > 0)
      {
	// Copy the face to a previous face
	surfmesh->face[n-num_removed].a = surfmesh->face[n].a;
	surfmesh->face[n-num_removed].b = surfmesh->face[n].b;
	surfmesh->face[n-num_removed].c = surfmesh->face[n].c;
	surfmesh->face[n-num_removed].m = surfmesh->face[n].m;
	surfmesh->face[n-num_removed].sel = surfmesh->face[n].sel;
      }

      // Un mark vertex for deletion
      surfmesh->vertex[surfmesh->face[n].a].m = 0;
      surfmesh->vertex[surfmesh->face[n].b].m = 0;
      surfmesh->vertex[surfmesh->face[n].c].m = 0;
    }
  }

  // Update the number of faces
  surfmesh->num_faces -= num_removed;
  SurfaceMesh_removeUnconnectedVertices(surfmesh);

}



/*
 * ***************************************************************************
 * SubRoutine: SurfaceMesh_getCenterRadius
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Return the center and radius of SurfaceMesh
 * ***************************************************************************
 */
ATOM SurfaceMesh_getCenterRadius(SurfaceMesh* surfmesh)
{
  ATOM data;
  int i;

  // Get the center and radius of the molecular surface mesh
  float mol_center_x = 0;
  float mol_center_y = 0;
  float mol_center_z = 0;
  float distance;
  
  for (i=0; i<surfmesh->num_vertices; i++)
  {
    mol_center_x += surfmesh->vertex[i].x;
    mol_center_y += surfmesh->vertex[i].y;
    mol_center_z += surfmesh->vertex[i].z;
  }
  if (surfmesh->num_vertices > 0) 
  {
    data.x = (float)(mol_center_x/(double)surfmesh->num_vertices);
    data.y = (float)(mol_center_y/(double)surfmesh->num_vertices);
    data.z = (float)(mol_center_z/(double)surfmesh->num_vertices);
  }
  else 
  {
    printf("no data found ...\n");
    data.radius = 0.0;
    data.x = 0.0;
    data.y = 0.0;
    data.z = 0.0;
    return data;
  }
  
  data.radius = 0;
  for (i=0; i < surfmesh->num_vertices; i++) 
  {
    distance = sqrt((surfmesh->vertex[i].x-data.x)*(surfmesh->vertex[i].x-data.x)+
		    (surfmesh->vertex[i].y-data.y)*(surfmesh->vertex[i].y-data.y)+
		    (surfmesh->vertex[i].z-data.z)*(surfmesh->vertex[i].z-data.z));
    if (distance > data.radius)
      data.radius = distance;
  }
  
  return data;
}

void SurfaceMesh_translate(SurfaceMesh* surfmesh, float dx, float dy, float dz)
{
  int i;
  for (i=0; i<surfmesh->num_vertices; i++)
  {
    surfmesh->vertex[i].x += dx;
    surfmesh->vertex[i].y += dy;
    surfmesh->vertex[i].z += dz;
  }
}

void SurfaceMesh_scale(SurfaceMesh* surfmesh, float scale_x, 
		       float scale_y, float scale_z)
{
  int i;
  for (i=0; i<surfmesh->num_vertices; i++)
  {
    surfmesh->vertex[i].x *= scale_x;
    surfmesh->vertex[i].y *= scale_y;
    surfmesh->vertex[i].z *= scale_z;
  }
}

void SurfaceMesh_scale(SurfaceMesh* surfmesh, float scale)
{
  SurfaceMesh_scale(surfmesh, scale, scale, scale);
}

void SurfaceMesh_centeralize(SurfaceMesh* surfmesh)
{
  ATOM center = SurfaceMesh_getCenterRadius(surfmesh);
  SurfaceMesh_translate(surfmesh, -center.x, -center.y, -center.z);
}

void SurfaceMesh_eigenvalues(SurfaceMesh* surfmesh){
  EIGENVECT eigen_vect;
  FLTVECT eigen_value;
  float max_angle;
  int n;

  // Check if neighborlist is created
  if (!surfmesh->neighbor_list)
    SurfaceMesh_createNeighborlist(surfmesh);

  if (surfmesh->neighbor_list == NULL)
  {
    printf("Could not create neighbor list some polygons might not be closed.\n");
    printf("Bailing out...\n");
    return ;
  }

  for (n = 0; n < surfmesh->num_vertices; n++)
  {

    // If we have a vertex wich is not selected we continue
    if (!surfmesh->vertex[n].sel)
      continue;
    
    eigen_vect = GetEigenVector(surfmesh, n, &eigen_value, &max_angle);

    // FIXME: Fix a way to get this information out in some way...
    printf("%d: (%7.2f, %7.2f, %7.2f), %5.2f, %5.2f, %5.2f\n", n, surfmesh->vertex[n].x,
	   surfmesh->vertex[n].y, surfmesh->vertex[n].z, eigen_value.x, eigen_value.y, eigen_value.z);
  }
  
}

void SurfaceMesh_splitMultipleConnectedSurfaces(SurfaceMesh* surfmesh)
{
  int m, n, a0, b0;
  int a, b, c, d;
  NPNT3 *first_ngr, *second_ngr, *tmp_ngr, *last_ngr;
  NPNT3 **neighbor_list;
  bool closed;

  // First stl library class... Let there be more...
  std::vector<int> non_manifold_vertices(0);

  // Re-create the neigborlist
  SurfaceMesh_createNeighborlist(surfmesh);
  
  // Nothing to do
  if (surfmesh->closed)
    return;

  // Get the neighbor list
  neighbor_list = surfmesh->neighbor_list;

  // Find non_manifold_vertices
  for (n = 0; n < surfmesh->num_vertices; n++) 
  {

    // Check that the neighbor list is connected
    tmp_ngr = neighbor_list[n];

    closed = true;
    while (tmp_ngr->next != NULL)
    {
      // Check that we are connected
      if (tmp_ngr->b != tmp_ngr->next->a)
      {
	if (std::find(non_manifold_vertices.begin(), non_manifold_vertices.end(), n) ==
	    non_manifold_vertices.end())
	{
	  non_manifold_vertices.push_back(n);
	}
	  
	// Do not bail, just register the vertex to not be done anything with
	surfmesh->vertex[n].sel = false;

	closed = false;
      }
      
      // Step one face forward
      tmp_ngr = tmp_ngr->next;
    }

    // Check if the list forms a closed ring
    if (closed && b0 != c)
    {
      closed = false;

      // FIXME: Should we add non_manifold vertex for this case too?
    }
    
  }

  int new_number_of_vertices=surfmesh->num_vertices;
  
  printf("Experimental method for splitting multiple connected surfaces mesh.\n");
  printf("Old number of vertices: %d\n", surfmesh->num_vertices);
  
  // Scratch data
  std::vector<NPNT3 *> tmp_neighborgs(0);
  std::vector<FLTVECT> tmp_vertices(0);
  std::vector<int> not_fixed_verts(0);

  printf("Trying to repair miss connected mesh.\n");
  printf("Number of non manifold vertices: %d,\n", (int)non_manifold_vertices.size());
  for (int j = 0; j < non_manifold_vertices.size(); j++)
  {
    // Get faulty vertex
    n = non_manifold_vertices[j];
    tmp_ngr = neighbor_list[n];
    a = tmp_ngr->a;
    bool fixed = true;
    int num_first = 0;
    int num_second = 0;

    // Gather data for averaging vert coordinates of connected vertices
    FLTVECT first_vert, second_vert;
    first_vert.x = first_vert.y = first_vert.z = 0;
    second_vert.x = second_vert.y = second_vert.z = 0;
    
    while (tmp_ngr->next != NULL)
    {
 	num_first++;

 	// Gather geometry information
 	first_vert.x += surfmesh->vertex[tmp_ngr->a].x;
 	first_vert.y += surfmesh->vertex[tmp_ngr->a].y;
 	first_vert.z += surfmesh->vertex[tmp_ngr->a].z;
 	
 	// If we have closed a ring but are not finished
 	if (tmp_ngr->b != tmp_ngr->next->a && tmp_ngr->b == a)
 	{
 	  
 	  //printf("Got one closed ring of %d faces for vert %d.\n", num_first, n);

 	  // Start a new neighbor ring
 	  last_ngr = tmp_ngr;
 	  tmp_ngr = tmp_ngr->next;
 	  
 	  // Push back the start of the next ring
 	  tmp_neighborgs.push_back(tmp_ngr);

 	  c = tmp_ngr->a;

 	  // Walk the next ring of connected vertices
 	  while(tmp_ngr->next != NULL)
 	  {
 	    num_second++;

 	    // Gather geometry information
 	    second_vert.x += surfmesh->vertex[tmp_ngr->a].x;
 	    second_vert.y += surfmesh->vertex[tmp_ngr->a].y;
 	    second_vert.z += surfmesh->vertex[tmp_ngr->a].z;

 	    // The ring has a hole
 	    if (tmp_ngr->b != tmp_ngr->next->a)
 	    {
 	      // If not already registered
 	      //if (std::find(not_fixed_verts.begin(), not_fixed_verts.end(), n) ==
 	      //	  not_fixed_verts.end())
 	      //	printf("Could not fix vertex %d.\n", n);
 	      //else
 	      not_fixed_verts.push_back(n);
 	      
 	      fixed = false;
 	      tmp_neighborgs[j] = NULL;
 	      break;
 	    }

 	    tmp_ngr = tmp_ngr->next;
 	    
 	  }
 	  
 	  // Add last neighbor
 	  num_second++;

 	  // Gather geometry information
 	  second_vert.x += surfmesh->vertex[tmp_ngr->a].x;
 	  second_vert.y += surfmesh->vertex[tmp_ngr->a].y;
 	  second_vert.z += surfmesh->vertex[tmp_ngr->a].z;

 	  // The ring is not closed
 	  if (tmp_ngr->b != c)
 	  {
 	    // If not already registered
 	    if (std::find(not_fixed_verts.begin(), not_fixed_verts.end(), n) ==
 		not_fixed_verts.end())
 	      printf("Could not fix vertex %d.\n", n);
 	    else
 	      not_fixed_verts.push_back(n);
 	    fixed = false;
 	    printf("Could not fix vertex %d.\n", n);
 	    tmp_neighborgs[j] = NULL;
 	  }
 	  
 	  // If we manage to fix it
 	  if (fixed)
 	  {
 	    //printf("Got a second closed ring of %d faces for vert %d.\n", num_second, n);

 	    // Bump the total number of vertices
 	    new_number_of_vertices += 1;

 	    // Select vertex
 	    surfmesh->vertex[n].sel = true;

 	    // Register the vertices for addition
 	    tmp_vertices.push_back(surfmesh->vertex[n]);
 	    FLTVECT& tmp_vertex = tmp_vertices[tmp_vertices.size()-1];

 	    // Find mid point of vertex rings
 	    first_vert.x /= num_first;
 	    first_vert.y /= num_first;
 	    first_vert.z /= num_first;
 	    second_vert.x /= num_second;
 	    second_vert.y /= num_second;
 	    second_vert.z /= num_second;

 	    // Find midpoint between original vertex and vertex ring and
 	    // update the two vertices with that point
 	    //printf("old coordinate: (%.2f, %.2f, %.2f)\n", surfmesh->vertex[n].x, 
 	    //	   surfmesh->vertex[n].y, surfmesh->vertex[n].z);
 	    
 	    surfmesh->vertex[n].x = (surfmesh->vertex[n].x+first_vert.x)/2;
 	    surfmesh->vertex[n].y = (surfmesh->vertex[n].y+first_vert.y)/2;
 	    surfmesh->vertex[n].z = (surfmesh->vertex[n].z+first_vert.z)/2;
 	    tmp_vertex.x = (tmp_vertex.x+second_vert.x)/2;
 	    tmp_vertex.y = (tmp_vertex.y+second_vert.y)/2;
 	    tmp_vertex.z = (tmp_vertex.z+second_vert.z)/2;
 	    

 	    //printf("New coordinate: (%.2f, %.2f, %.2f)\n", surfmesh->vertex[n].x, 
 	    //   surfmesh->vertex[n].y, surfmesh->vertex[n].z);
 	    //printf("New coordinate2: (%.2f, %.2f, %.2f)\n", tmp_vertex.x, 
 	    //	   tmp_vertex.y, tmp_vertex.z);

 	    // Update neighborlist of neighbors
 	    tmp_ngr = tmp_neighborgs[j];
 	    while (tmp_ngr != NULL)
 	    {
 	      // Update face information
 	      if (surfmesh->face[tmp_ngr->c].a == n)
 	      {
 	      	//printf("Changing connection for face %d.a (%d->%d)\n", tmp_ngr->c, 
 	      	//       n, new_number_of_vertices-1);
 	      	surfmesh->face[tmp_ngr->c].a = new_number_of_vertices-1;
 	      }
 	      
 	      if (surfmesh->face[tmp_ngr->c].b == n)
 	      {
 	      	//printf("Changing connection for face %d.b (%d->%d)\n", tmp_ngr->c, 
 	      	//       n, new_number_of_vertices-1);
 	      	surfmesh->face[tmp_ngr->c].b = new_number_of_vertices-1;
 	      }
 	      
 	      if (surfmesh->face[tmp_ngr->c].c == n)
 	      {
 	      	//printf("Changing connection for face %d.c (%d->%d)\n", tmp_ngr->c, 
 	      	//       n, new_number_of_vertices-1);
 	      	surfmesh->face[tmp_ngr->c].c = new_number_of_vertices-1;
 	      }

 	      tmp_ngr = tmp_ngr->next;
 	    }

 	    break;
 	  }
 	}
    
 	// Step one face forward
 	tmp_ngr = tmp_ngr->next;
 	
 	if (tmp_ngr == NULL)
 	  break;
    }
    
  }

  // If we have not fixed all vertices we check if some vertices share the 
  // same coordinate
  if (not_fixed_verts.size() > 0)
  {
    // Collect verts with same coordinate.
    std::vector<std::pair<int,int> > same_coord(0);
    float tol = 1e-3;
    for (int i=0; i<not_fixed_verts.size()-1; i++)
    {
 	for (int j=i+1; j<not_fixed_verts.size(); j++)
 	{
 	  if (std::fabs(surfmesh->vertex[not_fixed_verts[i]].x-
 			surfmesh->vertex[not_fixed_verts[j]].x)<tol&&
 	      std::fabs(surfmesh->vertex[not_fixed_verts[i]].y-
 			surfmesh->vertex[not_fixed_verts[j]].y)<tol&&
 	      std::fabs(surfmesh->vertex[not_fixed_verts[i]].z-
 			surfmesh->vertex[not_fixed_verts[j]].z)<tol)
 	  {
 	    same_coord.push_back(std::make_pair(not_fixed_verts[i], 
 						not_fixed_verts[j]));
 	    printf("Same coordinates: (%d, %d)\n", not_fixed_verts[i], 
 		   not_fixed_verts[j]);
 	  }
 	}
    }
    
  }
  
  printf("New number of vertices: %d\n", new_number_of_vertices);

  // No vertices added
  if (new_number_of_vertices == surfmesh->num_vertices)
    return;
  
  // Add collected vertices and neighbors
  FLTVECT* new_vertices = (FLTVECT*)malloc(sizeof(FLTVECT)*new_number_of_vertices);
  
  // Counter for added new entities
  m = 0;
  for (n = 0; n < new_number_of_vertices; n++)
  {
    
    // If we are not adding any vertices
    if (n<surfmesh->num_vertices)
    {
 	// Copy the old vertices
 	new_vertices[n].x = surfmesh->vertex[n].x;
 	new_vertices[n].y = surfmesh->vertex[n].y;
 	new_vertices[n].z = surfmesh->vertex[n].z;
 	new_vertices[n].sel = surfmesh->vertex[n].sel;
 	new_vertices[n].m = surfmesh->vertex[n].m;
 	
    }
    else
    {
 	
 	// Add new vertices
 	new_vertices[n].x = tmp_vertices[m].x;
 	new_vertices[n].y = tmp_vertices[m].y;
 	new_vertices[n].z = tmp_vertices[m].z;
 	new_vertices[n].sel = tmp_vertices[m].sel;
 	new_vertices[n].m = tmp_vertices[m].m;
 	
 	m++;
 	
    }
  }

  // Delete old data and reassign with new
  surfmesh->neighbor_list = neighbor_list;
  SurfaceMesh_destroyNeighborlist(surfmesh);

  free(surfmesh->vertex);
  surfmesh->vertex = new_vertices;
  surfmesh->num_vertices = new_number_of_vertices;
  surfmesh->closed = true;

  // Debug
  //SurfaceMesh_writeOFF(surfmesh, "new_surfmesh.off");
  
  // Re-create neighbors
  SurfaceMesh_createNeighborlist(surfmesh);

}

void SurfaceMesh_removeUnconnectedPatches(SurfaceMesh* surfmesh, int minimal_number)
{
  // Re-create the neigborlist
  SurfaceMesh_createNeighborlist(surfmesh);
  
  // Scratch data
  std::vector<int> will_visit(0);
  std::vector<int> patches_faces(surfmesh->num_faces, -1);
  std::vector<int> patches_vertices(surfmesh->num_vertices, -1);
  std::vector<bool> visited_verts(surfmesh->num_vertices, false);
  std::vector<bool> visited_faces(surfmesh->num_faces, false);
  will_visit.reserve(surfmesh->num_vertices);

  int patch = 0;
  int num_visited = 0;
  
  // Iterate over the different patches
  while (num_visited < surfmesh->num_vertices)
  {
    
    // Find next vertex which will be visited
    for (int i=0; i < surfmesh->num_vertices; i++)
    {

      // Check if we have visited the vert
      if (visited_verts[i])
	continue;
      
      // Schedule the vert to be visited
      will_visit.push_back(i);
      break;
    }
    
    // Iterate over all neigboring vertices in a patch
    while (will_visit.size()>0)
    {

      //  Bump number of visited verts
      num_visited += 1;

      // Vertex which we now visit
      const int visiting = will_visit.back();
      will_visit.pop_back();

      //printf("Patch: %d, Num visited: %d, Now visiting: %d, # to visit: %d\n", 
      //       patch, num_visited, visiting, will_visit.size());
      
      // Store the patch
      patches_vertices[visiting] = patch;

      // Get all neigbors
      NPNT3 *tmp_ngr = surfmesh->neighbor_list[visiting];
      while (tmp_ngr != NULL)
      {
	
	if (!visited_verts[tmp_ngr->a])
	{
	  will_visit.push_back(tmp_ngr->a);
	  visited_verts[tmp_ngr->a] = true;
	}
      
	patches_faces[tmp_ngr->c] = patch;
	tmp_ngr = tmp_ngr->next;
      }
    
    }

    // Bump the patch number
    patch += 1;
  }
  
  printf("Num patches: %d\n", patch);

  // Collect the number faces in each patch
  std::vector<int> num_faces_per_patch(patch, 0);
  for (int i=0; i<surfmesh->num_faces; i++)
  {
    num_faces_per_patch[patches_faces[i]] += 1;
  }

  for (int i = 0; i<patch; i++)
    printf("Num faces in patch %d: %d\n", i, num_faces_per_patch[i]);

}


