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


/* ***************************************************************************
 * File:     SurfSmooth.C   
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Smooth and coarsen surface triangular meshes
 * ****************************************************************************
 */

#include <gamer/biom.h>
#include "gamercf.h"

// Declare internal GAMer methods
FLTVECT GetPositionSurfaceOnly(float, float, float, int, int, int, SurfaceMesh *);
float GetAngleSurfaceOnly(SurfaceMesh *, int, int, int);
void  MoveVerticesSurfaceOnly(SurfaceMesh *, int);
FLTVECT GetNormals(SurfaceMesh *, int);
EIGENVECT GetEigenVector(SurfaceMesh *, int, FLTVECT *, float *);
void PolygonSubdivision(SurfaceMesh *, NPNT3 *, int *, int *, int);
char CheckFlipAction(SurfaceMesh *, int, int, int, int, bool);
void EdgeFlipping(SurfaceMesh *, int, bool);
float GetDotProduct(SurfaceMesh *, int, int, int);
FLTVECT Rotate(float, float, float, float, float, float);
FLTVECT GetCrossProduct(SurfaceMesh *, int, int, int);
void NormalSmooth(SurfaceMesh *, int);
FLTVECT GetNormals(SurfaceMesh *, int, int, int);

/*
 * ***************************************************************************
 * Routine:  GenerateHistogram   
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Generate the angle distribution (0 - 180 degrees)
 * ***************************************************************************
 */
void GenerateHistogram(SurfaceMesh *surfmesh)
{
  int n, m;
  int a, b, c;
  float angle;
  int histogram[18];

  for (m = 0; m < 18; m++)
    histogram[m] = 0;

  for (n = 0; n < surfmesh->num_faces; n++) 
  {
    a = surfmesh->face[n].a;
    b = surfmesh->face[n].b;
    c = surfmesh->face[n].c;

    angle = GetAngleSurfaceOnly(surfmesh, a, b, c);
    for (m=0; m<18; m++) 
    {
      if (angle >= m*10 && angle < m*10+10)
	histogram[m] += 1;
    }

    angle = GetAngleSurfaceOnly(surfmesh, b, a, c);
    for (m=0; m<18; m++) 
    {
      if (angle >= m*10 && angle < m*10+10)
        histogram[m] += 1;
    }

    angle = GetAngleSurfaceOnly(surfmesh, c, a, b);
    for (m=0; m<18; m++) 
    {
      if (angle >= m*10 && angle < m*10+10)
        histogram[m] += 1;
    }
   
  }
  
  for (m = 0; m < 18; m++) 
  {
    printf("%f  ", 100.0*(float)histogram[m]/((float)surfmesh->num_faces*3.0));
  }
  printf("\n\n");
}

/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_smooth   
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Function for surface smoothing
 * ***************************************************************************
 */
bool SurfaceMesh_smooth(SurfaceMesh *surfmesh, 
		unsigned int max_min_angle, 
		unsigned int min_max_angle, 
		unsigned int max_iter, 
		bool preserve_ridges)
{
  float min_angle, max_angle;
  int num_small, num_large;
  unsigned int i, n;
  bool smoothed = false;

  // Check if neighborlist is created
  if (surfmesh->neighbor_list == NULL)
    SurfaceMesh_createNeighborlist(surfmesh);
  
  if (surfmesh->neighbor_list == NULL)
  {
    printf("Could not create neigbor list some polygons might not be closed \nor you need to harmanize the normals of the mesh.\n");
    printf("Bailing out...\n");
    return false;
  }
  
  i = 0;

  SurfaceMesh_getMinMaxAngles(surfmesh, &min_angle, &max_angle, &num_small, 
			      &num_large, max_min_angle, min_max_angle);

  // Print the initial quality only when doing 1 or more iterations
  if (max_iter > 1) 
  {
    printf("Min Max angles\n");
    printf("%2d: min_angle: %f - max_angle: %f - "
	   "smaller-than-%d: %d - larger-than-%d: %d\n", 
	   i, min_angle, max_angle, max_min_angle, num_small, 
	   min_max_angle, num_large);

  }

  // Check if the mesh is smoothed
  smoothed = min_angle > max_min_angle && max_angle < min_max_angle;
  while (!smoothed && i < max_iter)
  {
    
    i++;
    
    // Smooth all vertices
    for (n = 0; n < surfmesh->num_vertices; n++)
    {

      // If we have a vertex wich is not selected we continue
      if (!surfmesh->vertex[n].sel)
	continue;

      MoveVerticesSurfaceOnly(surfmesh, n);
      EdgeFlipping(surfmesh, n, preserve_ridges);
    }

    // Calculate and print quality after surface smooth
    SurfaceMesh_getMinMaxAngles(surfmesh, &min_angle, &max_angle, &num_small, 
				&num_large, max_min_angle, min_max_angle);
    
    // Print the iteration number only when doing 1 or more iterations
    if (max_iter != 1)
      printf("%2d: min_angle: %f - max_angle: %f - "
	     "smaller-than-%d: %d - larger-than-%d: %d\n", 
	     i, min_angle, max_angle, max_min_angle, num_small, 
	     min_max_angle, num_large);
    else
      printf("    min_angle: %f - max_angle: %f - "
	     "smaller-than-%d: %d - larger-than-%d: %d\n", 
	     min_angle, max_angle, max_min_angle, num_small, 
	     min_max_angle, num_large);
      
    // Check if the mesh is smoothed
    smoothed = min_angle > max_min_angle && max_angle < min_max_angle;
  }

  return smoothed;
}

/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_normalSmooth
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Function for surface smoothing using normal smoothing
 * ***************************************************************************
 */
void SurfaceMesh_normalSmooth(SurfaceMesh *surfmesh)
{
  unsigned int n;
  float min_angle, max_angle;
  int num_small, num_large;
  
  // Check if neighborlist is created
  if (!surfmesh->neighbor_list)
    SurfaceMesh_createNeighborlist(surfmesh);

  if (surfmesh->neighbor_list == NULL)
  {
    printf("Could not create neigbor list some polygons might not be closed.\n");
    printf("Bailing out...\n");
    return ;
  }
    
  // Normal smooth all vertices
  for (n = 0; n < surfmesh->num_vertices; n++) 
  {
    
    // If we have a vertex wich is not selected we continue
    if (!surfmesh->vertex[n].sel)
    {
      continue;
    }

    NormalSmooth(surfmesh, n);
  }

  SurfaceMesh_getMinMaxAngles(surfmesh, &min_angle, &max_angle, 
			      &num_small, &num_large, 15, 150);
  printf("    min_angle: %f - max_angle: %f - smaller-than-15: %d "
	 "- larger-than-150: %d\n", 
	 min_angle, max_angle, num_small, num_large);
  
}


/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_assignActiveSites
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com), Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Assign markers to the vertices from a list of user-defined spheres.
 *           Note that vertices might be marked twice. Such vertex will then get
 *           the marker from the last sphere.
 * ***************************************************************************
 */
void SurfaceMesh_assignActiveSites(SurfaceMesh* surfmesh, ATOM* sphere_list, 
				   unsigned int num_spheres, 
				   unsigned int* sphere_markers)
{
  int i, n;
  float x, y, z;
  float center_x, center_y, center_z;
  float dist, radius;

  // Walk through all spheres and mark the vertices
  for (i = 0; i < num_spheres; i++) 
  {
    center_x = sphere_list[i].x;
    center_y = sphere_list[i].y;
    center_z = sphere_list[i].z;
    radius   = sphere_list[i].radius;

    for (n = 0; n < surfmesh->num_vertices; n++) 
    {
      x = surfmesh->vertex[n].x;
      y = surfmesh->vertex[n].y;
      z = surfmesh->vertex[n].z;
      dist = sqrt((x-center_x)*(x-center_x)+
		  (y-center_y)*(y-center_y)+
		  (z-center_z)*(z-center_z));

      // If the vertex n, is within the sphere mark it and deselect it 
      // for coarsening and smoothing
      if (dist < radius)
      {
	surfmesh->vertex[n].m = sphere_markers[i];
	surfmesh->vertex[n].sel = false;
	//printf("vertex: %d marked with %d\n", n, sphere_markers[i]);
      }
    }
  }
}

/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_coarse   
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Coarsen a surface mesh
 * ***************************************************************************
 */
char SurfaceMesh_coarse(SurfaceMesh* surfmesh, 
			float coarse_rate, 
			float flatness_rate, float denseness_weight, 
			float max_normal_angle)
{
  int m, n, a0, b0;
  int a, b, c, face_marker;
  float x, y, z;
  NPNT3 *first_ngr, *second_ngr, *tmp_ngr, *tmp_ngr1;
  int number, neighbor_number, num;
  float nx, ny, nz;
  EIGENVECT eigen_vect;
  FLTVECT eigen_value;
  float average_len, max_len;
  float ratio1 = 1.0, ratio2 = 1.0;
  int face_available_list[64], face_available_index;
  int neighbor_tmp_list[64];
  float weight, angle;

  int start_index, *vertex_index, *face_index;
 
  FLTVECT pos_vect;
  float w1, w2, w3;
  char delete_flag;
  float max_angle;
  
  char stop;
  int vertex_num;
  bool delete_vertex;

  // Check if neighborlist is created if not create it and reset vertex markers 
  // if they do not excist
  if (!surfmesh->neighbor_list)
    SurfaceMesh_createNeighborlist(surfmesh);

  if (surfmesh->neighbor_list == NULL)
  {
    printf("Could not create neigbor list some polygons might not be closed.\n");
    printf("Bailing out...\n");
    return 0;
  }
    
  NPNT3** neighbor_list = surfmesh->neighbor_list;

  vertex_index = (int *)malloc(sizeof(int)*surfmesh->num_vertices);
  face_index = (int *)malloc(sizeof(int)*surfmesh->num_faces);

  printf("Num vertices: %d, num faces: %d\n", surfmesh->num_vertices, 
	 surfmesh->num_faces);

  vertex_num = surfmesh->num_vertices;

  stop = 0;
  // If using sparseness weight, calculate the average segment length of the mesh
  if (denseness_weight > 0.0)
  {
    average_len = 0;
    for (n = 0; n < surfmesh->num_faces; n++) 
    {
      a = surfmesh->face[n].a;
      b = surfmesh->face[n].b;
      c = surfmesh->face[n].c;
     
      nx = sqrt((surfmesh->vertex[a].x-surfmesh->vertex[b].x)*
		(surfmesh->vertex[a].x-surfmesh->vertex[b].x)+
		(surfmesh->vertex[a].y-surfmesh->vertex[b].y)*
		(surfmesh->vertex[a].y-surfmesh->vertex[b].y)+
		(surfmesh->vertex[a].z-surfmesh->vertex[b].z)*
		(surfmesh->vertex[a].z-surfmesh->vertex[b].z));
      ny = sqrt((surfmesh->vertex[a].x-surfmesh->vertex[c].x)*
		(surfmesh->vertex[a].x-surfmesh->vertex[c].x)+
		(surfmesh->vertex[a].y-surfmesh->vertex[c].y)*
		(surfmesh->vertex[a].y-surfmesh->vertex[c].y)+
		(surfmesh->vertex[a].z-surfmesh->vertex[c].z)*
		(surfmesh->vertex[a].z-surfmesh->vertex[c].z));
      nz = sqrt((surfmesh->vertex[c].x-surfmesh->vertex[b].x)*
		(surfmesh->vertex[c].x-surfmesh->vertex[b].x)+
		(surfmesh->vertex[c].y-surfmesh->vertex[b].y)*
		(surfmesh->vertex[c].y-surfmesh->vertex[b].y)+
		(surfmesh->vertex[c].z-surfmesh->vertex[b].z)*
		(surfmesh->vertex[c].z-surfmesh->vertex[b].z));
      average_len += (nx+ny+nz)/3.0f;
    }
    if (surfmesh->num_faces == 0) 
    {
      printf("Zero degree on a vertex ....\n");
      printf("Bailing out...\n");
      return 0;
    }
    else 
      surfmesh->avglen = average_len/(float)(surfmesh->num_faces);
  }

  // The main loop over all vertices
  for (n = 0; n < surfmesh->num_vertices; n++) 
  {

    // Status report
    if (((n+1) % 888) == 0 || (n+1) == surfmesh->num_vertices) 
    {
      printf("%2.2f%% done (%08d)          \r", 100.0*(n+1)/
	     (float)surfmesh->num_vertices, n+1);
      fflush(stdout);
    }

    // If the vertex have been flagged to not be removed
    if (!surfmesh->vertex[n].sel)
    {
      //printf("Do not remove vertex %d\n", n);
      continue;
    }

    // Check if the vertex has enough neigborgs to be deleted
    delete_flag = 1;
    first_ngr = neighbor_list[n];
    fflush(stdout);
    while (first_ngr != NULL) 
    {
      a = first_ngr->a;
      number = 0;
      num = 0;
      second_ngr = neighbor_list[a];
      while (second_ngr != NULL) 
      {
	fflush(stdout);
	b = second_ngr->a;
	tmp_ngr = neighbor_list[n];
        while (tmp_ngr != NULL) {
          if (tmp_ngr->a == b)
            num++;
          tmp_ngr = tmp_ngr->next;
        }
	number++;
	second_ngr = second_ngr->next;
      }
      
      if (number <= 3 || num > 2)
	delete_flag = 0;
      first_ngr = first_ngr->next;
    }
    
    if (delete_flag)
    {
      x = surfmesh->vertex[n].x;
      y = surfmesh->vertex[n].y;
      z = surfmesh->vertex[n].z;
      
      max_len = -1;
      first_ngr = neighbor_list[n];

      // If using sparseness as a criteria for coarsening
      // calculate the maximal segment length
      if (denseness_weight > 0.0)
      {
        while (first_ngr != NULL) 
	{
	  a = first_ngr->a;
	  b = first_ngr->b;
	  
	  nx = sqrt((x-surfmesh->vertex[a].x)*(x-surfmesh->vertex[a].x)+
		    (y-surfmesh->vertex[a].y)*(y-surfmesh->vertex[a].y)+
		    (z-surfmesh->vertex[a].z)*(z-surfmesh->vertex[a].z));
	  ny = sqrt((x-surfmesh->vertex[b].x)*(x-surfmesh->vertex[b].x)+
		    (y-surfmesh->vertex[b].y)*(y-surfmesh->vertex[b].y)+
		    (z-surfmesh->vertex[b].z)*(z-surfmesh->vertex[b].z));
	  if (nx > max_len)
	    max_len = nx;
	  if (ny > max_len)
	    max_len = ny;
	  
	  first_ngr = first_ngr->next;
        }

	// Max segment length over the average segment length of the mesh
	ratio2 = max_len/surfmesh->avglen;
	ratio2 = pow(ratio2, denseness_weight);
      }
      
      // If using curvatory as a coarsening criteria
      // calculate the local structure tensor
      if (flatness_rate > 0.0) 
      {
        eigen_vect = GetEigenVector(surfmesh, n, &eigen_value, &max_angle);
	
	if (eigen_value.x == 0) 
	{
	  printf("max eigen_value is zero.... \n");
	  exit(0);
	}
	else 
	{
	  //printf("Eigenvalues: %d, (%.3f, %.3f, %.3f)\n", n,
	  //	 eigen_value.x, eigen_value.y, eigen_value.z);
	  ratio1 = fabs((eigen_value.y)/(eigen_value.x));
	  ratio1 = pow(ratio1, flatness_rate);
	  //ratio1 = (1.0-max_angle)*fabs((eigen_value.y)/(eigen_value.x));
	}
      }
      
      // Compare the two coarseness criterias against the given coarse_rate
      delete_vertex = ratio1*ratio2 < coarse_rate;

      // Use maximal angle between vertex normal as a complementary coarse criteria
      if (max_normal_angle > 0)
	delete_vertex = delete_vertex && max_angle > max_normal_angle;
      
      // Deleting a vertex and retrianglulate the hole
      if (delete_vertex) 
      {
	vertex_num--;
	
	/* delete vertex n */
	surfmesh->vertex[n].x = -99999;
	surfmesh->vertex[n].y = -99999;
	surfmesh->vertex[n].z = -99999;
	
	neighbor_number = 0;
	first_ngr = neighbor_list[n];
	while (first_ngr != NULL) 
	{
	  a = first_ngr->a;
	  c = first_ngr->c;
	  face_available_list[neighbor_number] = c;
	  neighbor_tmp_list[neighbor_number] = a;
	  neighbor_number++;
	  
	  // Get face marker
	  face_marker = surfmesh->face[c].m;

	  /* delete faces associated with vertex n */
	  surfmesh->face[c].a = -1;
	  surfmesh->face[c].b = -1;
	  surfmesh->face[c].c = -1;
	  surfmesh->face[c].m = -1;
	  
	  /* delete neighbors associated with vertex n */
	  second_ngr = neighbor_list[a];
	  tmp_ngr = second_ngr;
	  while (second_ngr != NULL) 
	  {
	    if (second_ngr->a == n || second_ngr->b == n) 
	    {
	      if (second_ngr == neighbor_list[a]) 
	      {
		neighbor_list[a] = second_ngr->next;
		free(second_ngr);
		second_ngr = neighbor_list[a];
		tmp_ngr = second_ngr;
	      }
	      else 
	      {
		tmp_ngr->next = second_ngr->next;
		free(second_ngr);
     		second_ngr = tmp_ngr->next;
	      }
	    }
	    else 
	    {
	      if (second_ngr == neighbor_list[a]) 
	      {
		second_ngr = second_ngr->next;
	      }
	      else 
	      {
		tmp_ngr = second_ngr;
		second_ngr = second_ngr->next;
	      }
	    }
	  }
	  
	  number = 0;
	  second_ngr = neighbor_list[a];
	  while (second_ngr != NULL) 
	  {
	    number++;
	    second_ngr = second_ngr->next;
	  }
	  first_ngr->b = number;
	  first_ngr = first_ngr->next;
	}
	
	first_ngr = neighbor_list[n];
	while (first_ngr->next != NULL) 
	  first_ngr = first_ngr->next;
	first_ngr->next = neighbor_list[n];
	
	face_available_index = 0;
	PolygonSubdivision(surfmesh, neighbor_list[n], face_available_list, 
			   &face_available_index, face_marker);
	
	/* order the neighbors */
	for (m = 0; m < neighbor_number; m++) 
	{
	  first_ngr = neighbor_list[neighbor_tmp_list[m]];
	  c = first_ngr->a;
	  while (first_ngr != NULL) 
	  {
	    a = first_ngr->a;
	    b = first_ngr->b;

	    second_ngr = first_ngr->next;
	    while (second_ngr != NULL) 
	    {
	      a0 = second_ngr->a;
	      b0 = second_ngr->b;

	      // Assume counter clockwise orientation
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
	    if (first_ngr->next == NULL) 
	    {
	      if (first_ngr->b != c) 
	      {
		printf("some polygons are not closed: %d \n", n);
		// exit(0);
	      }
	    }
	    
	    first_ngr = first_ngr->next;
	  }
	}
	
	/* Smooth the neighbors */
	for (m = 0; m < neighbor_number; m++) 
	{
	  if (!surfmesh->vertex[num].sel)
	    continue;
	  num = neighbor_tmp_list[m];
	  x = surfmesh->vertex[num].x;
	  y = surfmesh->vertex[num].y;
	  z = surfmesh->vertex[num].z;
	  nx = 0;
	  ny = 0;
	  nz = 0;
	  
	  weight = 0;
	  first_ngr = neighbor_list[num];
	  while (first_ngr != NULL) 
	  {
	    a = first_ngr->a;
	    b = first_ngr->b;
	    second_ngr = first_ngr->next;
	    if (second_ngr == NULL)
	      second_ngr = neighbor_list[num];
	    c = second_ngr->b;
	    pos_vect = GetPositionSurfaceOnly(x, y, z, b, a, c, surfmesh);
	    angle = GetDotProduct(surfmesh, b, a, c);
	    angle += 1.0;
	    nx += angle*pos_vect.x;
	    ny += angle*pos_vect.y;
	    nz += angle*pos_vect.z;

	    weight += angle;
	    first_ngr = first_ngr->next;
	  }
	  
	  if (weight > 0) 
	  {
	    nx /= weight;
	    ny /= weight;
	    nz /= weight;

	    eigen_vect = GetEigenVector(surfmesh, num, &eigen_value, &max_angle);
	    if ((eigen_vect.x1==0 && eigen_vect.y1==0 && eigen_vect.z1==0) ||
		(eigen_vect.x2==0 && eigen_vect.y2==0 && eigen_vect.z2==0) ||
		(eigen_vect.x3==0 && eigen_vect.y3==0 && eigen_vect.z3==0)) 
	    {
	      surfmesh->vertex[num].x = nx;
	      surfmesh->vertex[num].y = ny;
	      surfmesh->vertex[num].z = nz;
	    }
	    else 
	    {
	      nx -= x;
	      ny -= y;
	      nz -= z;
	      w1 = (nx*eigen_vect.x1+ny*eigen_vect.y1+nz*eigen_vect.z1)/
		(1.0+eigen_value.x);
	      w2 = (nx*eigen_vect.x2+ny*eigen_vect.y2+nz*eigen_vect.z2)/
		(1.0+eigen_value.y);
	      w3 = (nx*eigen_vect.x3+ny*eigen_vect.y3+nz*eigen_vect.z3)/
		(1.0+eigen_value.z);
	      surfmesh->vertex[num].x = w1*eigen_vect.x1+w2*eigen_vect.x2+
		w3*eigen_vect.x3 + x;
	      surfmesh->vertex[num].y = w1*eigen_vect.y1+w2*eigen_vect.y2+
		w3*eigen_vect.y3 + y;
	      surfmesh->vertex[num].z = w1*eigen_vect.z1+w2*eigen_vect.z2+
		w3*eigen_vect.z3 + z;
	    }
	  }
	}
      }
    }
    /*
      if (vertex_num < MeshSizeUpperLimit) {
      stop = 1;
      break;
      }
    */
  }
  
  /* Clean the lists of nodes and faces */
  start_index = 0;
  for (n = 0; n < surfmesh->num_vertices; n++) 
  {
    if (surfmesh->vertex[n].x != -99999 &&
	surfmesh->vertex[n].y != -99999 &&
	surfmesh->vertex[n].z != -99999) 
    {
      if (start_index != n) 
      {
	surfmesh->vertex[start_index].x = surfmesh->vertex[n].x;
	surfmesh->vertex[start_index].y = surfmesh->vertex[n].y;
	surfmesh->vertex[start_index].z = surfmesh->vertex[n].z;
	surfmesh->vertex[start_index].m = surfmesh->vertex[n].m;
	surfmesh->vertex[start_index].sel = surfmesh->vertex[n].sel;
	neighbor_list[start_index] = neighbor_list[n];
      }
      
      vertex_index[n] = start_index;
      start_index++;
    }
    else 
    {
      vertex_index[n] = -1;
    }
  }

  surfmesh->num_vertices = start_index;

  start_index = 0;
  for (n = 0; n < surfmesh->num_faces; n++) 
  {
    a = surfmesh->face[n].a;
    b = surfmesh->face[n].b;
    c = surfmesh->face[n].c;
    face_marker = surfmesh->face[n].m;
    if (a >= 0 && b >= 0 && c >= 0 && 
	vertex_index[a] >= 0 && vertex_index[b] >= 0 && vertex_index[c] >= 0) 
    {
      surfmesh->face[start_index].a = vertex_index[a];
      surfmesh->face[start_index].b = vertex_index[b];
      surfmesh->face[start_index].c = vertex_index[c];
      surfmesh->face[start_index].m = face_marker;
      face_index[n] = start_index;
      start_index++;
    }
    else 
    {
      face_index[n] = -1;
    }
  }
  surfmesh->num_faces = start_index;
  
  for (n = 0; n < surfmesh->num_vertices; n++) 
  {
    first_ngr = neighbor_list[n];
    while (first_ngr != NULL) 
    {
      a = first_ngr->a;
      b = first_ngr->b;
      c = first_ngr->c;
    
      first_ngr->a = vertex_index[a];
      first_ngr->b = vertex_index[b];
      first_ngr->c = face_index[c];
      
      first_ngr = first_ngr->next;
    }
  }
  
  free(vertex_index);
  free(face_index);
  printf("\n");
  return(stop);
}

/*
 * ***************************************************************************
 * Routine:  PolygonSubdivision
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Recursively re-triangulate the "empty" polygon
 * ***************************************************************************
 */
void PolygonSubdivision(SurfaceMesh *surfmesh, 
			NPNT3 *start_ngr, int *face_available_list, 
			int *face_available_index, int face_marker)
{
  NPNT3 **neighbor_list = surfmesh->neighbor_list;
  NPNT3 *first_ngr, *second_ngr;
  NPNT3 *tmp_ngr, *first_copy_ngr, *second_copy_ngr;
  int min_num, degree;
  int face_index, number;
  int a, b, c;

  
  number = 1;
  tmp_ngr = start_ngr;
  while (tmp_ngr->next != start_ngr) 
  {
    number++;
    tmp_ngr = tmp_ngr->next;
  }
   
  if (number < 3) 
  {
    printf("error: number of nodes less than 3 \n");
    exit(0);
  }

  if (number == 3) 
  {
    a = start_ngr->a;
    tmp_ngr = start_ngr->next;
    free(start_ngr);
    start_ngr = tmp_ngr;

    b = start_ngr->a;
    tmp_ngr = start_ngr->next;
    free(start_ngr);
    start_ngr = tmp_ngr;

    c = start_ngr->a;
    tmp_ngr = start_ngr->next;
    free(start_ngr);
    start_ngr = tmp_ngr;

    face_index = face_available_list[*face_available_index];
    surfmesh->face[face_index].a = a;
    surfmesh->face[face_index].b = b;
    surfmesh->face[face_index].c = c;
    surfmesh->face[face_index].m = face_marker;
    *face_available_index += 1;

    first_ngr = (NPNT3 *)malloc(sizeof(NPNT3));
    first_ngr->a = b;
    first_ngr->b = c;
    first_ngr->c = face_index;
    first_ngr->next = neighbor_list[a];
    neighbor_list[a] = first_ngr;

    first_ngr = (NPNT3 *)malloc(sizeof(NPNT3));
    first_ngr->a = c;
    first_ngr->b = a;
    first_ngr->c = face_index;
    first_ngr->next = neighbor_list[b];
    neighbor_list[b] = first_ngr;

    first_ngr = (NPNT3 *)malloc(sizeof(NPNT3));
    first_ngr->a = a;
    first_ngr->b = b;
    first_ngr->c = face_index;
    first_ngr->next = neighbor_list[c];
    neighbor_list[c] = first_ngr;

  }
  else 
  {
    tmp_ngr = start_ngr;
    min_num = tmp_ngr->b;
    first_ngr = tmp_ngr;
    tmp_ngr = tmp_ngr->next;
    while (tmp_ngr != start_ngr) 
    {
      degree = tmp_ngr->b;
      if (degree < min_num) 
      {
	min_num = degree;
	first_ngr = tmp_ngr;
      }
      tmp_ngr = tmp_ngr->next;
    }
    
    min_num = 99999;
    tmp_ngr = start_ngr;
    if (tmp_ngr != first_ngr &&
	tmp_ngr != first_ngr->next &&
	tmp_ngr->next != first_ngr) 
    {
      min_num = tmp_ngr->b;
      second_ngr = tmp_ngr;
    }
    
    tmp_ngr = tmp_ngr->next;
    while (tmp_ngr != start_ngr) 
    {
      degree = tmp_ngr->b;
      if (tmp_ngr != first_ngr &&
	  tmp_ngr != first_ngr->next &&
	  tmp_ngr->next != first_ngr &&
	  degree < min_num) 
      {
	min_num = degree;
	second_ngr = tmp_ngr;
      }
      tmp_ngr = tmp_ngr->next;
    }

    first_ngr->b += 1;
    second_ngr->b += 1;
    first_copy_ngr = (NPNT3 *)malloc(sizeof(NPNT3));
    first_copy_ngr->a = first_ngr->a;
    first_copy_ngr->b = first_ngr->b;
    second_copy_ngr = (NPNT3 *)malloc(sizeof(NPNT3));
    second_copy_ngr->a = second_ngr->a;
    second_copy_ngr->b = second_ngr->b;
    tmp_ngr = first_ngr;
    while (tmp_ngr->next != first_ngr)
      tmp_ngr = tmp_ngr->next;
    tmp_ngr->next = first_copy_ngr;
    first_copy_ngr->next = second_copy_ngr;
    second_copy_ngr->next = second_ngr->next;
    second_ngr->next = first_ngr;

    PolygonSubdivision(surfmesh, first_ngr, face_available_list, 
		       face_available_index, face_marker);
    PolygonSubdivision(surfmesh, first_copy_ngr, face_available_list, 
		       face_available_index, face_marker);
  }

  return;
}




/*
 * ***************************************************************************
 * Routine:  MoveVerticesSurfaceOnly
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Smooth the surface mesh by moving each of the vertices
 *           using a combination of angle-based method 
 *           and local structure tensor 
 * ***************************************************************************
 */
void MoveVerticesSurfaceOnly(SurfaceMesh *surfmesh, int n)
{
  int a, b, c;
  float x, y, z;
  FLTVECT pos_vect;
  NPNT3 **neighbor_list = surfmesh->neighbor_list;
  NPNT3 *first_ngr, *second_ngr;
  float weight, angle;
  float nx, ny, nz;
  EIGENVECT eigen_vect;
  FLTVECT eigen_value;
  float w1, w2, w3;
  float max_angle;
  
  
  x = surfmesh->vertex[n].x;
  y = surfmesh->vertex[n].y;
  z = surfmesh->vertex[n].z;
  
  nx = 0;
  ny = 0;
  nz = 0;
  
  weight = 0;
  first_ngr = neighbor_list[n];
  while (first_ngr != NULL) 
  {
    a = first_ngr->a;
    b = first_ngr->b;
    second_ngr = first_ngr->next;
    if (second_ngr == NULL)
      second_ngr = neighbor_list[n];
    c = second_ngr->b;
    pos_vect = GetPositionSurfaceOnly(x, y, z, b, a, c, surfmesh);
    angle = GetDotProduct(surfmesh, b, a, c);
    angle += 1.0;
    nx += angle*pos_vect.x;
    ny += angle*pos_vect.y;
    nz += angle*pos_vect.z;
    
    weight += angle;
    first_ngr = first_ngr->next;
  }

  if (weight > 0) 
  {
    nx /= weight;
    ny /= weight;
    nz /= weight;

    eigen_vect = GetEigenVector(surfmesh, n, &eigen_value, &max_angle);
    if ((eigen_vect.x1==0 && eigen_vect.y1==0 && eigen_vect.z1==0) ||
	(eigen_vect.x2==0 && eigen_vect.y2==0 && eigen_vect.z2==0) ||
	(eigen_vect.x3==0 && eigen_vect.y3==0 && eigen_vect.z3==0)) 
    {
      //printf("old point (%0.2f, %0.2f, %0.2f), new point (%0.2f, %0.2f, %0.2f)\n", 
      //	     surfmesh->vertex[n].x, surfmesh->vertex[n].y, surfmesh->vertex[n].z, 
      //	     nx, ny, nz);
      surfmesh->vertex[n].x = nx;
      surfmesh->vertex[n].y = ny;
      surfmesh->vertex[n].z = nz;
      
    }
    else 
    {
      nx -= x;
      ny -= y;
      nz -= z;
      w1 = (nx*eigen_vect.x1+ny*eigen_vect.y1+nz*eigen_vect.z1)/(1.0+eigen_value.x);
      w2 = (nx*eigen_vect.x2+ny*eigen_vect.y2+nz*eigen_vect.z2)/(1.0+eigen_value.y);
      w3 = (nx*eigen_vect.x3+ny*eigen_vect.y3+nz*eigen_vect.z3)/(1.0+eigen_value.z);
      nx = w1*eigen_vect.x1+w2*eigen_vect.x2+w3*eigen_vect.x3 + x;
      ny = w1*eigen_vect.y1+w2*eigen_vect.y2+w3*eigen_vect.y3 + y;
      nz = w1*eigen_vect.z1+w2*eigen_vect.z2+w3*eigen_vect.z3 + z;

      //printf("old point (%0.2f, %0.2f, %0.2f), new point (%0.2f, %0.2f, %0.2f)\n", 
      //	     surfmesh->vertex[n].x, surfmesh->vertex[n].y, surfmesh->vertex[n].z, 
      //	     nx, ny, nz);

      surfmesh->vertex[n].x = nx;
      surfmesh->vertex[n].y = ny;
      surfmesh->vertex[n].z = nz;
    }
  }
}

/*
 * ***************************************************************************
 * Routine:  GetPositionSurfaceOnly   
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Move each of the vertices based on the angle-based method 
 * ***************************************************************************
 */
FLTVECT GetPositionSurfaceOnly(float x, float y, float z, int a, int b, int c, SurfaceMesh *surfmesh)
{
  float ax, ay, az;
  float bx, by, bz;
  float cx, cy, cz;
  float xx, yy, zz;
  float distance;
  FLTVECT tmp;
  float tx, ty, tz;

  
  ax = surfmesh->vertex[a].x;
  ay = surfmesh->vertex[a].y;
  az = surfmesh->vertex[a].z;
  bx = surfmesh->vertex[b].x;
  by = surfmesh->vertex[b].y;
  bz = surfmesh->vertex[b].z;
  cx = surfmesh->vertex[c].x;
  cy = surfmesh->vertex[c].y;
  cz = surfmesh->vertex[c].z;
  
  bx -= ax;
  by -= ay;
  bz -= az;
  distance = sqrt(bx*bx+by*by+bz*bz);
  if (distance > 0) 
  {
    bx /= distance;
    by /= distance;
    bz /= distance;
  }
  cx -= ax;
  cy -= ay;
  cz -= az;
  distance = sqrt(cx*cx+cy*cy+cz*cz);
  if (distance > 0) 
  {
    cx /= distance;
    cy /= distance;
    cz /= distance;
  }
  tx = 0.5*(cx+bx);
  ty = 0.5*(cy+by);
  tz = 0.5*(cz+bz);
  distance = sqrt(tx*tx+ty*ty+tz*tz);
  if (distance > 0) 
  {
    tx /= distance;
    ty /= distance;
    tz /= distance;
  }
  xx = by*cz-bz*cy;
  yy = bz*cx-bx*cz;
  zz = bx*cy-by*cx;
  distance = sqrt(xx*xx+yy*yy+zz*zz);
  if (distance > 0) 
  {
    xx /= distance;
    yy /= distance;
    zz /= distance;
  }
  bx = xx;
  by = yy;
  bz = zz;
  
  distance = tx*(x-ax)+ty*(y-ay)+tz*(z-az);
  xx = distance*tx + ax;
  yy = distance*ty + ay;
  zz = distance*tz + az;

  distance = bx*(x-xx)+by*(y-yy)+bz*(z-zz);
  tmp.x = distance*bx + xx;
  tmp.y = distance*by + yy;
  tmp.z = distance*bz + zz;

  return(tmp);
}
  

/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_getMinMaxAngles
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Calculate the minimum and maximum angles in a surface mesh 
 * ***************************************************************************
 */
void SurfaceMesh_getMinMaxAngles(SurfaceMesh *surfmesh, float *minangle, 
				 float *maxangle, int *num_small, int *num_large, 
				 int max_min_angle, int min_max_angle)
{
  int n, num1, num2;
  int a, b, c;
  float min_angle, max_angle;
  float angle;


  min_angle = 99999.0;
  max_angle = -99999.0;
  num1 = 0;
  num2 = 0;
  for (n = 0; n < surfmesh->num_faces; n++)
  {
    a = surfmesh->face[n].a;
    b = surfmesh->face[n].b;
    c = surfmesh->face[n].c;
     
    angle = GetAngleSurfaceOnly(surfmesh, a, b, c);
    if (angle != -999) 
    {
      if (angle < min_angle)
	min_angle = angle;
      if (angle > max_angle)
	max_angle = angle;
      if (angle < max_min_angle)
	num1++;
      if (angle > min_max_angle)
	num2++;
    }
    angle = GetAngleSurfaceOnly(surfmesh, b, a, c);
    if (angle != -999) 
    {
      if (angle < min_angle)
	min_angle = angle;
      if (angle > max_angle)
	max_angle = angle;
      if (angle < max_min_angle)
        num1++;
      if (angle > min_max_angle)
        num2++;
    }
    angle = GetAngleSurfaceOnly(surfmesh, c, a, b);
    if (angle != -999) 
    {
      if (angle < min_angle)
	min_angle = angle;
      if (angle > max_angle)
	max_angle = angle;
      if (angle < max_min_angle)
        num1++;
      if (angle > min_max_angle)
        num2++;
    }
  }
  
  *minangle = min_angle;
  *maxangle = max_angle;
  *num_small = num1;
  *num_large = num2;
}

/*
 * ***************************************************************************
 * Routine:  GetAngleSurfaceOnly   
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Calculate the angle defined by three vertices 
 * ***************************************************************************
 */
float GetAngleSurfaceOnly(SurfaceMesh *surfmesh, int a, int b, int c)
{
  float ax, ay, az;
  float bx, by, bz;
  float cx, cy, cz;
  float length1, length2, length3;
  float angle;
  
  ax = surfmesh->vertex[a].x;
  ay = surfmesh->vertex[a].y;
  az = surfmesh->vertex[a].z;
  bx = surfmesh->vertex[b].x;
  by = surfmesh->vertex[b].y;
  bz = surfmesh->vertex[b].z;
  cx = surfmesh->vertex[c].x;
  cy = surfmesh->vertex[c].y;
  cz = surfmesh->vertex[c].z;

  length1 = (ax-bx)*(ax-bx)+(ay-by)*(ay-by)+(az-bz)*(az-bz);
  length2 = (ax-cx)*(ax-cx)+(ay-cy)*(ay-cy)+(az-cz)*(az-cz);
  length3 = (bx-cx)*(bx-cx)+(by-cy)*(by-cy)+(bz-cz)*(bz-cz);
  if (length1 == 0 || length2 == 0)
    angle = -999;
  else 
  {
    angle = 0.5*(length1+length2-length3)/sqrt(length1*length2);
    angle = acos(angle)*180.0/PIE;
  }

  return(angle);
}


/*
 * ***************************************************************************
 * Routine:  GetEigenVector
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Calculate the eigenvalues and eigenvectors of 
 *           the local structure tensor in a neighborhood of "index0" 
 * ***************************************************************************
 */
EIGENVECT GetEigenVector(SurfaceMesh *surfmesh, 
			 int index0, FLTVECT *eigen_value, float *max_ang)
{
  int index, dist;
  int n, m;
  double x1, x2, x3;
  double a, b, Q;
  double c0, c1, c2;
  double A[3][3];
  double B[6];
  double theta, p;
  double tx, ty, tz;
  EIGENVECT tmp;
  FLTVECT normal, normal0;
  
  int IndexArray[333];
  int DistArray[333];
  int start_ptr, end_ptr;
  int visited;

  NPNT3 **neighbor_list = surfmesh->neighbor_list;
  NPNT3 *first_ngr;
  float angle, max_angle;
  
  normal = GetNormals(surfmesh, index0);
  
  //printf("Normal: %d (%.2f, %.2f, %.2f)\n", index0, normal.x, normal.y, normal.z);
  A[0][0] = normal.x*normal.x;
  A[0][1] = normal.x*normal.y;
  A[0][2] = normal.x*normal.z;
  A[1][1] = normal.y*normal.y;
  A[1][2] = normal.y*normal.z;
  A[2][2] = normal.z*normal.z;

  start_ptr = 0;
  end_ptr = 1;
  IndexArray[start_ptr] = index0;
  DistArray[start_ptr] = 0;

  max_angle = 99999.0;
  normal0.x = normal.x;
  normal0.y = normal.y;
  normal0.z = normal.z;
  while (start_ptr < end_ptr) 
  {
    index = IndexArray[start_ptr];
    dist = DistArray[start_ptr];
    start_ptr ++;

    if (dist < ((DIM_SCALE>2) ? (3):(2))) 
    {
      first_ngr = neighbor_list[index];
      while (first_ngr != NULL) 
      {
	m = first_ngr->a;
	visited = 0;
	for (n = 0; n < end_ptr; n++) 
	{
	  if (IndexArray[n] == m) 
	  {
	    visited = 1;
	    //printf("Visited!\n");
	    break;
	  }
	}
	if (visited == 0) 
	{
	  normal = GetNormals(surfmesh, m);
	  angle = normal0.x*normal.x+normal0.y*normal.y+normal0.z*normal.z;
	  if (angle < 0)
	    angle = -angle;
	  if (angle < max_angle)
	    max_angle = angle;
	  A[0][0] += normal.x*normal.x;
	  A[0][1] += normal.x*normal.y;
	  A[0][2] += normal.x*normal.z;
	  A[1][1] += normal.y*normal.y;
	  A[1][2] += normal.y*normal.z;
	  A[2][2] += normal.z*normal.z;
	  IndexArray[end_ptr] = m;
	  DistArray[end_ptr] = dist+1;
	  end_ptr ++;
	}
	first_ngr = first_ngr->next;
      }
    }
  }
  *max_ang = max_angle;

  A[1][0] = A[0][1];
  A[2][0] = A[0][2];
  A[2][1] = A[1][2];
  
  c0 = A[0][0]*A[1][1]*A[2][2]+2*A[0][1]*A[0][2]*A[1][2]-A[0][0]*A[1][2]*A[1][2]
    -A[1][1]*A[0][2]*A[0][2]-A[2][2]*A[0][1]*A[0][1];
  c1 = A[0][0]*A[1][1]-A[0][1]*A[0][1]+A[0][0]*A[2][2]-
    A[0][2]*A[0][2]+A[1][1]*A[2][2]-A[1][2]*A[1][2];
  c2 = A[0][0]+A[1][1]+A[2][2];
  
  a = (3.0*c1-c2*c2)/3.0;
  b = (-2.0*c2*c2*c2+9.0*c1*c2-27.0*c0)/27.0;
  Q = b*b/4.0+a*a*a/27.0;
  	    
  theta = atan2(sqrt(-Q), -0.5*b);
  p = sqrt(0.25*b*b-Q);
  
  x1 = c2/3.0+2.0*pow(p, 1.0/3.0)*cos(theta/3.0);
  x2 = c2/3.0-pow(p, 1.0/3.0)*(cos(theta/3.0)+sqrt(3.0)*sin(theta/3.0));
  x3 = c2/3.0-pow(p, 1.0/3.0)*(cos(theta/3.0)-sqrt(3.0)*sin(theta/3.0));

  // If we have a perfect flat area inline of one of the x, y, z axis
  if (isnan(x1) || isnan(x2) || isnan(x3)) 
  {
    
    //printf("Got a NAN: %d\n", index0);
    eigen_value->x = c2;
    eigen_value->y = 0; 
    eigen_value->z = 0; 
    
    tmp.x1 = 1;
    tmp.y1 = 0;
    tmp.z1 = 0;

    tmp.x2 = 0;
    tmp.y2 = 1;
    tmp.z2 = 0;

    tmp.x3 = 0;
    tmp.y3 = 0;
    tmp.z3 = 1;

    return tmp;
  }
  
  tx = max(x1, max(x2, x3));
  if (tx == x1) 
  {
    if (x2 >= x3) 
    {
  	ty = x2;
  	tz = x3;
    }
    else 
    {
  	ty = x3;
  	tz = x2;
    }
  }
  else if (tx == x2) 
  {
    if (x1 >= x3) 
    {
  	ty = x1;
  	tz = x3;
    }
    else {
  	ty = x3;
  	tz = x1;
    }
  }
  else if (tx == x3) 
  {
    if (x1 >= x2) 
    {
  	ty = x1;
  	tz = x2;
    }
    else 
    {
  	ty = x2;
  	tz = x1;
    }
  }
  x1 = tx;
  x2 = ty;
  x3 = tz;
  eigen_value->x = tx;
  eigen_value->y = ty;
  eigen_value->z = tz;
  
  if (x1 > 99999 || x1 < -99999 ||
  	x2 > 99999 || x2 < -99999 ||
  	x3 > 99999 || x3 < -99999) 
  {
    printf("dsadsadsad: %f %f %f\n", x1, x2, x3);
    exit(0);
  }

  
  A[0][0] -= x1;
  A[1][1] -= x1;
  A[2][2] -= x1;
  B[0] = A[1][1]*A[2][2]-A[1][2]*A[1][2];
  B[1] = A[0][2]*A[1][2]-A[0][1]*A[2][2];
  B[2] = A[0][0]*A[2][2]-A[0][2]*A[0][2];
  B[3] = A[0][1]*A[1][2]-A[0][2]*A[1][1];
  B[4] = A[0][1]*A[0][2]-A[1][2]*A[0][0];
  B[5] = A[0][0]*A[1][1]-A[0][1]*A[0][1];
  c0 = B[0]*B[0]+B[1]*B[1]+B[3]*B[3];
  c1 = B[1]*B[1]+B[2]*B[2]+B[4]*B[4];
  c2 = B[3]*B[3]+B[4]*B[4]+B[5]*B[5];
  if (c0 >= c1 && c0 >= c2) 
  {
    tx = B[0];
    ty = B[1];
    tz = B[3];
  }
  else if (c1 >= c0 && c1 >= c2) 
  {
    tx = B[1];
    ty = B[2];
    tz = B[4];
  }
  else if (c2 >= c0 && c2 >= c1) 
  {
    tx = B[3];
    ty = B[4];
    tz = B[5];
  }
  p = sqrt(tx*tx+ty*ty+tz*tz);
  if (p > 0) 
  {
    tx /= p;
    ty /= p;
    tz /= p;
  }
  tmp.x1 = tx;
  tmp.y1 = ty;
  tmp.z1 = tz;
  A[0][0] += x1;
  A[1][1] += x1;
  A[2][2] += x1;
  
  
  A[0][0] -= x2;
  A[1][1] -= x2;
  A[2][2] -= x2;
  B[0] = A[1][1]*A[2][2]-A[1][2]*A[1][2];
  B[1] = A[0][2]*A[1][2]-A[0][1]*A[2][2];
  B[2] = A[0][0]*A[2][2]-A[0][2]*A[0][2];
  B[3] = A[0][1]*A[1][2]-A[0][2]*A[1][1];
  B[4] = A[0][1]*A[0][2]-A[1][2]*A[0][0];
  B[5] = A[0][0]*A[1][1]-A[0][1]*A[0][1];
  c0 = B[0]*B[0]+B[1]*B[1]+B[3]*B[3];
  c1 = B[1]*B[1]+B[2]*B[2]+B[4]*B[4];
  c2 = B[3]*B[3]+B[4]*B[4]+B[5]*B[5];
  if (c0 >= c1 && c0 >= c2) 
  {
    tx = B[0];
    ty = B[1];
    tz = B[3];
  }
  else if (c1 >= c0 && c1 >= c2) 
  {
    tx = B[1];
    ty = B[2];
    tz = B[4];
  }
  else if (c2 >= c0 && c2 >= c1) 
  {
    tx = B[3];
    ty = B[4];
    tz = B[5];
  }
  p = sqrt(tx*tx+ty*ty+tz*tz);
  if (p > 0) 
  {
    tx /= p;
    ty /= p;
    tz /= p;
  }
  tmp.x2 = tx;
  tmp.y2 = ty;
  tmp.z2 = tz;
  
  tmp.x3 = tmp.y1*tz-tmp.z1*ty;
  tmp.y3 = tmp.z1*tx-tmp.x1*tz;
  tmp.z3 = tmp.x1*ty-tmp.y1*tx;

  return(tmp);
}

/*
 * ***************************************************************************
 * Routine:  GetNormals   
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Calculate the normal vector at vertex "n" 
 * ***************************************************************************
 */
FLTVECT GetNormals(SurfaceMesh *surfmesh, int n)
{
  int a, b;
  float x, y, z;
  NPNT3 **neighbor_list = surfmesh->neighbor_list;
  NPNT3 *first_ngr;
  int number;
  FLTVECT normal;
  float gx, gy, gz;
  float ax, ay, az;
  float bx, by, bz;
  float length;
 
  x = surfmesh->vertex[n].x;
  y = surfmesh->vertex[n].y;
  z = surfmesh->vertex[n].z;
  
  number = 0;
  normal.x = 0;
  normal.y = 0;
  normal.z = 0;
  first_ngr = neighbor_list[n];
  while (first_ngr != NULL) 
  {
    a = first_ngr->a;
    b = first_ngr->b;
    
    ax = surfmesh->vertex[a].x-x;
    ay = surfmesh->vertex[a].y-y;
    az = surfmesh->vertex[a].z-z;
    length = sqrt(ax*ax+ay*ay+az*az);
    if (length > 0) 
    {
      ax /= length;
      ay /= length;
      az /= length;
    }
    bx = surfmesh->vertex[b].x-x;
    by = surfmesh->vertex[b].y-y;
    bz = surfmesh->vertex[b].z-z;
    length = sqrt(bx*bx+by*by+bz*bz);
    if (length > 0) 
    {
      bx /= length;
      by /= length;
      bz /= length;
    }
    gx = ay*bz-az*by;
    gy = az*bx-ax*bz;
    gz = ax*by-ay*bx;
    length = sqrt(gx*gx+gy*gy+gz*gz);
    if (length > 0) 
    {
      gx /= length;
      gy /= length;
      gz /= length;
    }
    length = normal.x*gx+normal.y*gy+normal.z*gz;
    if (length < 0) 
    {
      gx = -gx;
      gy = -gy;
      gz = -gz;
    }
    normal.x += gx;
    normal.y += gy;
    normal.z += gz;
    
    number ++;
    first_ngr = first_ngr->next;
  }
  
  if (number > 0) 
  {
    normal.x /= (float)number;
    normal.y /= (float)number;
    normal.z /= (float)number;
    length = sqrt(normal.x*normal.x+normal.y*normal.y+normal.z*normal.z);
    if (length > 0) 
    {
      normal.x /= length;
      normal.y /= length;
      normal.z /= length;
    }
  }
  else 
  {
    normal.x = 0;
    normal.y = 0;
    normal.z = 0;
  }

  return(normal);
}




/*
 * ***************************************************************************
 * Routine:  CheckFlipAction
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Check if the edge flipping is needed or not, by the 
 *           "smaller angle criterion" 
 * ***************************************************************************
 */
char CheckFlipAction(SurfaceMesh *surfmesh, 
		     int a, int b, int c, int d, bool preserve_ridges)
{
  NPNT3 **neighbor_list = surfmesh->neighbor_list;
  float min_angle1, min_angle2, angle;
  FLTVECT normal_a, normal_b;
  
  /* smaller angle criterion */
  min_angle1 = -99999;
  angle = GetDotProduct(surfmesh, a, b, c);
  if (angle > min_angle1)
    min_angle1 = angle;
  angle = GetDotProduct(surfmesh, a, b, d);
  if (angle > min_angle1)
    min_angle1 = angle;
  angle = GetDotProduct(surfmesh, b, a, c);
  if (angle > min_angle1)
    min_angle1 = angle;
  angle = GetDotProduct(surfmesh, b, a, d);
  if (angle > min_angle1)
    min_angle1 = angle;
  
  min_angle2 = -99999;
  angle = GetDotProduct(surfmesh, c, a, d);
  if (angle > min_angle2)
    min_angle2 = angle;
  angle = GetDotProduct(surfmesh, c, b, d);
  if (angle > min_angle2)
    min_angle2 = angle;
  angle = GetDotProduct(surfmesh, d, a, c);
  if (angle > min_angle2)
    min_angle2 = angle;
  angle = GetDotProduct(surfmesh, d, b, c);
  if (angle > min_angle2)
    min_angle2 = angle;
  
  // Check which of the triangle combination has the smallest angle
  // min_angle1 is the minimal angle of the flipped configuration
  // min_angle2 is the minimal angle of the present configuration
  if (min_angle1 > min_angle2)
  {
    // Check if the angle between the normals of the two triangles are 
    // too small for a flip action, for example if we are on a ridge
    normal_a = GetCrossProduct(surfmesh, a, c, b);
    normal_b = GetCrossProduct(surfmesh, a, b, d);
    
    // If we want to preserve the ridges the angle between the surface
    // normals must be  smaller than cos(60deg) to flip the edges
    if (not preserve_ridges || 
	normal_a.x*normal_b.x+normal_a.y*normal_b.y+normal_a.z*normal_b.z>0.866)
      return(1);
  }
  return(0);
}


/*
 * ***************************************************************************
 * Routine:  GetDotProduct   
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Calculate the dot product of two vectors (b-a) and (c-a) 
 * ***************************************************************************
 */
float GetDotProduct(SurfaceMesh *surfmesh, int a, int b, int c)
{
  float cx, cy, cz;
  float bx, by, bz;
  float length;


  bx = surfmesh->vertex[b].x-surfmesh->vertex[a].x;
  by = surfmesh->vertex[b].y-surfmesh->vertex[a].y;
  bz = surfmesh->vertex[b].z-surfmesh->vertex[a].z;
  length = sqrt(bx*bx+by*by+bz*bz);
  if (length > 0) 
  {
    bx /= length;
    by /= length;
    bz /= length;
  }
  
  cx = surfmesh->vertex[c].x-surfmesh->vertex[a].x;
  cy = surfmesh->vertex[c].y-surfmesh->vertex[a].y;
  cz = surfmesh->vertex[c].z-surfmesh->vertex[a].z;
  length = sqrt(cx*cx+cy*cy+cz*cz);
  if (length > 0) 
  {
    cx /= length;
    cy /= length;
    cz /= length;
  }
  
  length = bx*cx+by*cy+bz*cz;
  return(length);
}


/*
 * ***************************************************************************
 * Routine:  EdgeFlipping
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Perform the edge flipping 
 * ***************************************************************************
 */
void EdgeFlipping(SurfaceMesh *surfmesh, int n, bool preserve_ridges)
{
  int a, b, c;
  NPNT3 **neighbor_list = surfmesh->neighbor_list;
  NPNT3 *first_ngr, *second_ngr;
  NPNT3 *tmp_ngr1, *tmp_ngr2, *tmp_ngr;
  char flip_flag, flip_check;
  int f1, f2, number;
  float ax, ay, az;
  
  first_ngr = neighbor_list[n];
  while (first_ngr != NULL) 
  {

    number = 0;
    tmp_ngr = neighbor_list[n];
    while (tmp_ngr != NULL) 
    {
      number++;
      tmp_ngr = tmp_ngr->next;
    }
    if (number <= 3) 
    {
      if (number > 0) 
      {
	ax = 0;
	ay = 0;
	az = 0;
	tmp_ngr = neighbor_list[n];
	while (tmp_ngr != NULL) 
	{
	  a = tmp_ngr->a;
	  ax += surfmesh->vertex[a].x;
	  ay += surfmesh->vertex[a].y;
	  az += surfmesh->vertex[a].z;
	  tmp_ngr = tmp_ngr->next;
	}
	
	surfmesh->vertex[n].x = ax/(float)number;
	surfmesh->vertex[n].y = ay/(float)number;
	surfmesh->vertex[n].z = az/(float)number;
      }
      return;
    }

    a = first_ngr->a;
    b = first_ngr->b;
    second_ngr = first_ngr->next;
    if (second_ngr == NULL)
      second_ngr = neighbor_list[n];
    c = second_ngr->b;

    flip_flag = 1;
    number = 0;
    tmp_ngr = neighbor_list[b];
    while (tmp_ngr != NULL) 
    {
      number++;
      tmp_ngr = tmp_ngr->next;
    }
    if (number <= 3)
      flip_flag = 0;

    tmp_ngr = neighbor_list[a];
    while (tmp_ngr != NULL) 
    {
      if (tmp_ngr->a == c)
	flip_flag = 0;
      tmp_ngr = tmp_ngr->next;
    }
    tmp_ngr = neighbor_list[c];
    while (tmp_ngr != NULL) 
    {
      if (tmp_ngr->a == a)
	flip_flag = 0;
      tmp_ngr = tmp_ngr->next;
    }

    if (flip_flag) {
    
      flip_check = CheckFlipAction(surfmesh, n, b, a, c, preserve_ridges);
      
      if (flip_check) {
	f1 = first_ngr->c;
	f2 = second_ngr->c;
	
	/* Update face info */
	surfmesh->face[f1].a = n;
	surfmesh->face[f1].b = a; 
	surfmesh->face[f1].c = c; 
	surfmesh->face[f2].a = b; 
	surfmesh->face[f2].b = c; // Switch a and c here to make the face 
	surfmesh->face[f2].c = a; // normal point outward		  
	
	/* Delete the entries in neighbor lists */
	first_ngr->b = c;
	if (first_ngr->next == NULL) 
	  neighbor_list[n] = neighbor_list[n]->next;
	else
	  first_ngr->next = second_ngr->next;
	tmp_ngr1 = second_ngr;

	tmp_ngr = neighbor_list[b];
	while (tmp_ngr != NULL) 
	{
	  if (tmp_ngr->b == n) 
	    break;
	  tmp_ngr = tmp_ngr->next;
	}
	if (tmp_ngr == NULL)
	  printf("my god ... %d\n", n);
	if (tmp_ngr->a == c) 
	{
	  tmp_ngr->b = a;
	  tmp_ngr->c = f2;
	  if (tmp_ngr->next == NULL) 
	  {
	    second_ngr = neighbor_list[b];
	    neighbor_list[b] = second_ngr->next;
	  }
	  else 
	  {
	    second_ngr = tmp_ngr->next;
	    tmp_ngr->next = second_ngr->next;
	  }
	  tmp_ngr2 = second_ngr;
	}
	else 
	{
	  printf("delete error!!! %d : %d %d %d\n", n, a, b, c);
	  printf("(%f, %f, %f)\n", 
		 surfmesh->vertex[n].x, 
		 surfmesh->vertex[n].y, 
		 surfmesh->vertex[n].z);
	}
      
	/* Add the entries in neighbor lists */
	tmp_ngr = neighbor_list[a];
	while (tmp_ngr != NULL) 
	{
	  if ((tmp_ngr->a == n && tmp_ngr->b == b) ||
	      (tmp_ngr->a == b && tmp_ngr->b == n))
	    break;
	  tmp_ngr = tmp_ngr->next;
	}

	// Assume neigbors are stored counter clockwise
	if (tmp_ngr->a == b && tmp_ngr->b == n) 
	{
	  tmp_ngr->b = c;
	  tmp_ngr->c = f2;
	  tmp_ngr1->a = c;
	  tmp_ngr1->b = n;
	  tmp_ngr1->c = f1;
	  tmp_ngr1->next = tmp_ngr->next;
	  tmp_ngr->next = tmp_ngr1;
	}
	else 
	  printf("add error 111\n");
	
	tmp_ngr = neighbor_list[c];
	while (tmp_ngr != NULL)
	{
	  if ((tmp_ngr->a == n && tmp_ngr->b == b) ||
	      (tmp_ngr->a == b && tmp_ngr->b == n))
	    break;
	  tmp_ngr = tmp_ngr->next;
	}

	// Assume neigbors are stored counter clockwise
	if (tmp_ngr->a == n && tmp_ngr->b == b) 
	{
	  tmp_ngr->b = a;
	  tmp_ngr->c = f1;
	  tmp_ngr2->a = a;
	  tmp_ngr2->b = b;
	  tmp_ngr2->c = f2;
	  tmp_ngr2->next = tmp_ngr->next;
	  tmp_ngr->next = tmp_ngr2;
	}
	else 
	  printf("add error 222\n");
      }
    }
    
    first_ngr = first_ngr->next;
  }
  
}

/*
 * ***************************************************************************
 * Routine:  NormalSmooth
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Smooth the surface mesh by anisotropic normal-based averaging 
 * ***************************************************************************
 */
void NormalSmooth(SurfaceMesh *surfmesh, int n)
{
  int a, b, c, d, e;
  NPNT3 **neighbor_list = surfmesh->neighbor_list;
  NPNT3 *first_ngr, *second_ngr, *third_ngr;
  NPNT3 *tmp_ngr;
  float bx, by, bz;
  float cx, cy, cz;
  float dx, dy, dz;
  float fx, fy, fz;
  float gx, gy, gz;
  float pos_x, pos_y, pos_z;
  int number, num;
  float theta, phi, alpha;
  float length;
  FLTVECT normal, sv;
  
  
  number = 0;
  pos_x = 0;
  pos_y = 0;
  pos_z = 0;
  first_ngr = neighbor_list[n];
  while (first_ngr != NULL) 
  {
    a = first_ngr->a;
    b = first_ngr->b;
    second_ngr = first_ngr->next;
    if (second_ngr == NULL)
      second_ngr = neighbor_list[n];
    c = second_ngr->b;
    third_ngr = second_ngr->next;
    if (third_ngr == NULL)
      third_ngr = neighbor_list[n];
    d = third_ngr->b;
    
    tmp_ngr = neighbor_list[b];

    // If a vertex is neigbor with a non selected vertex continue
    if (!surfmesh->vertex[b].sel)
      return;
    
    while (tmp_ngr != NULL) 
    {
      if ((tmp_ngr->a == c && tmp_ngr->b != n) ||
	  (tmp_ngr->b == c && tmp_ngr->a != n))
	break;
      tmp_ngr = tmp_ngr->next;
    }
    if (tmp_ngr->a == c && tmp_ngr->b != n)
      e = tmp_ngr->b;
    else if (tmp_ngr->b == c && tmp_ngr->a != n)
      e = tmp_ngr->a;
    else 
      printf("normal smoothing error...\n");

    normal = GetCrossProduct(surfmesh, n, b, c);
    gx = normal.x;
    gy = normal.y;
    gz = normal.z;
    dx = 0;
    dy = 0;
    dz = 0;

    num  = 0;
    normal = GetCrossProduct(surfmesh, n, a, b);
    length = normal.x*gx+normal.y*gy+normal.z*gz;
    if (length > 0) 
    {
      num++;
      dx += length*normal.x;
      dy += length*normal.y;
      dz += length*normal.z;
    }
    normal = GetCrossProduct(surfmesh, n, c, d);
    length = normal.x*gx+normal.y*gy+normal.z*gz;
    if (length > 0) 
    {
      num++;
      dx += length*normal.x;
      dy += length*normal.y;
      dz += length*normal.z;
    }
    normal = GetCrossProduct(surfmesh, b, e, c);
    length = normal.x*gx+normal.y*gy+normal.z*gz;
    if (length > 0) 
    {
      num++;
      dx += length*normal.x;
      dy += length*normal.y;
      dz += length*normal.z;
    }

    length = sqrt(dx*dx+dy*dy+dz*dz);
    if (length > 0) 
    {
      dx /= length;
      dy /= length;
      dz /= length;
      fx = gy*dz-gz*dy;
      fy = gz*dx-gx*dz;
      fz = gx*dy-gy*dx;
      cx = surfmesh->vertex[c].x;
      cy = surfmesh->vertex[c].y;
      cz = surfmesh->vertex[c].z;
      bx = surfmesh->vertex[b].x;
      by = surfmesh->vertex[b].y;
      bz = surfmesh->vertex[b].z;
      length = fx*(bx-cx)+fy*(by-cy)+fz*(bz-cz);
      if (length >= 0) 
      {
	theta = (float)atan2(by-cy, bx-cx);
	phi = (float)atan2(bz-cz, sqrt((bx-cx)*(bx-cx)+(by-cy)*(by-cy)));
      }
      else 
      {
	theta = (float)atan2(cy-by, cx-bx);
	phi = (float)atan2(cz-bz, sqrt((bx-cx)*(bx-cx)+(by-cy)*(by-cy)));
      }
      
      alpha = acos(dx*gx+dy*gy+dz*gz)/(float)(4.0-num);
      sv = Rotate(surfmesh->vertex[n].x-cx, surfmesh->vertex[n].y-cy, surfmesh->vertex[n].z-cz, theta, phi, alpha);
      pos_x += sv.x+cx;
      pos_y += sv.y+cy;
      pos_z += sv.z+cz;
      
      number++;
    }
    
    first_ngr = first_ngr->next;
  }
  
  if (number > 0 && !isnan(pos_x) && !isnan(pos_y) && !isnan(pos_z)) 
  {
    surfmesh->vertex[n].x = pos_x/(float)number;
    surfmesh->vertex[n].y = pos_y/(float)number;
    surfmesh->vertex[n].z = pos_z/(float)number;
  }
  
}



/*
 * ***************************************************************************
 * Routine:  GetCrossProduct   
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Calculate the cross product vector between (b-a) and (c-a) 
 * ***************************************************************************
 */
FLTVECT GetCrossProduct(SurfaceMesh *surfmesh, int a, int b, int c) 
{
  float gx, gy, gz;
  float cx, cy, cz;
  float bx, by, bz;
  float length;
  FLTVECT value;


  bx = surfmesh->vertex[b].x-surfmesh->vertex[a].x;
  by = surfmesh->vertex[b].y-surfmesh->vertex[a].y;
  bz = surfmesh->vertex[b].z-surfmesh->vertex[a].z;
  length = sqrt(bx*bx+by*by+bz*bz);
  if (length > 0) 
  {
    bx /= length;
    by /= length;
    bz /= length;
  }
  cx = surfmesh->vertex[c].x-surfmesh->vertex[a].x;
  cy = surfmesh->vertex[c].y-surfmesh->vertex[a].y;
  cz = surfmesh->vertex[c].z-surfmesh->vertex[a].z;
  length = sqrt(cx*cx+cy*cy+cz*cz);
  if (length > 0) 
  {
    cx /= length;
    cy /= length;
    cz /= length;
  }
  gx = cy*bz-cz*by;
  gy = cz*bx-cx*bz;
  gz = cx*by-cy*bx;
  length = sqrt(gx*gx+gy*gy+gz*gz);
  if (length > 0) 
  {
    gx /= length;
    gy /= length;
    gz /= length;
  }

  value.x = gx;
  value.y = gy;
  value.z = gz;

  return(value);
}


/*
 * ***************************************************************************
 * Routine:  Rotate   
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Rotate a point "sx, sy, sz" around "theta, phi" axis by "angle" 
 * ***************************************************************************
 */
FLTVECT Rotate(float sx, float sy, float sz, 
	       float theta, float phi, float angle)
{
  float x, y, z;
  float xx, yy, zz;
  float a[3][3], b[3][3];
  FLTVECT tmp;


  a[0][0] = (float)(cos(0.5*PIE-phi)*cos(theta));
  a[0][1] = (float)(cos(0.5*PIE-phi)*sin(theta));
  a[0][2] = (float)-sin(0.5*PIE-phi);
  a[1][0] = (float)-sin(theta);
  a[1][1] = (float)cos(theta);
  a[1][2] = 0.f;
  a[2][0] = (float)(sin(0.5*PIE-phi)*cos(theta));
  a[2][1] = (float)(sin(0.5*PIE-phi)*sin(theta));
  a[2][2] = (float)cos(0.5*PIE-phi);

  b[0][0] = (float)(cos(0.5*PIE-phi)*cos(theta));
  b[0][1] = (float)-sin(theta); 
  b[0][2] = (float)(sin(0.5*PIE-phi)*cos(theta)); 
  b[1][0] = (float)(cos(0.5*PIE-phi)*sin(theta));
  b[1][1] = (float)cos(theta);
  b[1][2] = (float)(sin(0.5*PIE-phi)*sin(theta));
  b[2][0] = (float)-sin(0.5*PIE-phi);
  b[2][1] = 0.f;
  b[2][2] = (float)cos(0.5*PIE-phi);


  x = a[0][0]*sx+a[0][1]*sy+a[0][2]*sz;
  y = a[1][0]*sx+a[1][1]*sy+a[1][2]*sz;
  z = a[2][0]*sx+a[2][1]*sy+a[2][2]*sz;
      
  xx = (float)(cos(angle)*x - sin(angle)*y);
  yy = (float)(sin(angle)*x + cos(angle)*y);
  zz = z;

  tmp.x = b[0][0]*xx+b[0][1]*yy+b[0][2]*zz;
  tmp.y = b[1][0]*xx+b[1][1]*yy+b[1][2]*zz;
  tmp.z = b[2][0]*xx+b[2][1]*yy+b[2][2]*zz;
  
  return(tmp);
  
}
