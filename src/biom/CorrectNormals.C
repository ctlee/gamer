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

#include <gamer/biom.h>

/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_correctNormals   
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Function to correct normals so all facet normals point outwards
 * ***************************************************************************
 */
void SurfaceMesh_correctNormals(SurfaceMesh* surfmesh)
{
  int n, m;
  int i, j;
  int a, b, c, d;
  int a0, b0, c0;
  int *stack;
  unsigned char *visited;
  int start, end;
  float x, y, z;
  float distance, max_dist;
  int *face_list;
  int progress;
  float max_x, max_y, max_z;
  float ax, ay, az;
  float bx, by, bz;
  float cx, cy, cz;

  // get the max coordinates of the nodes
  max_x = -999999.0f;
  max_y = -999999.0f;
  max_z = -999999.0f;
  for (n = 0; n < surfmesh->num_vertices; n++) {
    x = surfmesh->vertex[n].x;
    y = surfmesh->vertex[n].y;
    z = surfmesh->vertex[n].z;
    if (x > max_x)
      max_x = x;
    if (y > max_y)
      max_y = y;
    if (z > max_z)
      max_z = z;
  }
  
  //printf("max coordinates: %f %f %f\n",max_x,max_y,max_z);
  max_x += 10.0f;
  max_y += 10.0f;
  max_z += 10.0f;

  // update the normal directions
  stack = (int*)malloc(sizeof(int)*surfmesh->num_faces);
  visited = (unsigned char*)malloc(sizeof(unsigned char)*surfmesh->num_faces);
  for (i = 0; i < surfmesh->num_faces; i++)
    visited[i] = 0;
  start = 0;
  progress = 0;
  while (1) {
    // find the unvisited face that is closest to (max_x,max_y,max_z)
    max_dist = 999999.0f;
    for (i = 0; i < surfmesh->num_faces; i++) {
      if (visited[i] == 0) {
	a = surfmesh->face[i].a;
	b = surfmesh->face[i].b;
	c = surfmesh->face[i].c;
	x = (surfmesh->vertex[a].x+surfmesh->vertex[b].x+surfmesh->vertex[c].x)/3.0f;
	y = (surfmesh->vertex[a].y+surfmesh->vertex[b].y+surfmesh->vertex[c].y)/3.0f;
	z = (surfmesh->vertex[a].z+surfmesh->vertex[b].z+surfmesh->vertex[c].z)/3.0f;
	distance = sqrt((x-max_x)*(x-max_x)+(y-max_y)*(y-max_y)+(z-max_z)*(z-max_z));
	if (distance < max_dist) {
	  j = i;
	  max_dist = distance;
	}
      }
    }
    if (max_dist == 999999.0f)
      break;
    
    a = surfmesh->face[j].a;
    b = surfmesh->face[j].b;
    c = surfmesh->face[j].c;
    ax = surfmesh->vertex[a].x;
    ay = surfmesh->vertex[a].y;
    az = surfmesh->vertex[a].z;
    bx = surfmesh->vertex[b].x;
    by = surfmesh->vertex[b].y;
    bz = surfmesh->vertex[b].z;
    cx = surfmesh->vertex[c].x;
    cy = surfmesh->vertex[c].y;
    cz = surfmesh->vertex[c].z;
    x = (ax+bx+cx)/3.0f;
    y = (ay+by+cy)/3.0f;
    z = (az+bz+cz)/3.0f;
    bx -= ax;
    by -= ay;
    bz -= az;
    cx -= ax;
    cy -= ay;
    cz -= az;
    ax = by*cz-bz*cy;
    ay = bz*cx-bx*cz;
    az = bx*cy-by*cx;
    if (ax*(max_x-x)+ay*(max_y-y)+az*(max_z-z) < 0) {
      // switch the order of the last two vertices
      b = surfmesh->face[j].b;
      surfmesh->face[j].b = surfmesh->face[j].c;
      surfmesh->face[j].c = b;
    }
    
    progress += start;
    stack[0] = j;
    start = 0;
    end = 1;  
    while (start < end) {
      if (((start) % 500) == 0) {
	printf("%2.2f%% done (%08d)\r", 100.0*(float)(start+progress)/(float)(surfmesh->num_faces), start+progress);
	fflush(stdout);
      }
      
      i = stack[start];
      visited[i] = 1;
      start++;
      a0 = surfmesh->face[i].a;
      b0 = surfmesh->face[i].b;
      c0 = surfmesh->face[i].c;
      
      for (n = 0; n < surfmesh->num_faces; n++) {
	if (visited[n] != 1) {
	  a = surfmesh->face[n].a;
	  b = surfmesh->face[n].b;
	  c = surfmesh->face[n].c;
	  
	  if (visited[n] == 0) {
	    if (a == a0 && b == b0) {
	      j = surfmesh->face[n].a;
	      surfmesh->face[n].a = surfmesh->face[n].b;
	      surfmesh->face[n].b = j;
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (b == a0 && a == b0) {
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (c == a0 && a == b0) {
	      j = surfmesh->face[n].a;
	      surfmesh->face[n].a = surfmesh->face[n].c;
	      surfmesh->face[n].c = j;
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (a == a0 && c == b0) {
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (b == a0 && c == b0) {
	      j = surfmesh->face[n].b;
	      surfmesh->face[n].b = surfmesh->face[n].c;
	      surfmesh->face[n].c = j;
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (c == a0 && b == b0) {
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    
	    else if (b == a0 && a == c0) {
	      j = surfmesh->face[n].a;
	      surfmesh->face[n].a = surfmesh->face[n].b;
	      surfmesh->face[n].b = j;
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (a == a0 && b == c0) {
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (a == a0 && c == c0) {
	      j = surfmesh->face[n].a;
	      surfmesh->face[n].a = surfmesh->face[n].c;
	      surfmesh->face[n].c = j;
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (c == a0 && a == c0) {
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (c == a0 && b == c0) {
	      j = surfmesh->face[n].b;
	      surfmesh->face[n].b = surfmesh->face[n].c;
	      surfmesh->face[n].c = j;
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (b == a0 && c == c0) {
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    
	    else if (a == b0 && b == c0) {
	      j = surfmesh->face[n].a;
	      surfmesh->face[n].a = surfmesh->face[n].b;
	      surfmesh->face[n].b = j;
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (b == b0 && a == c0) {
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (c == b0 && a == c0) {
	      j = surfmesh->face[n].a;
	      surfmesh->face[n].a = surfmesh->face[n].c;
	      surfmesh->face[n].c = j;
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (a == b0 && c == c0) {
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (b == b0 && c == c0) {
	      j = surfmesh->face[n].b;
	      surfmesh->face[n].b = surfmesh->face[n].c;
	      surfmesh->face[n].c = j;
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	    else if (c == b0 && b == c0) {
	      visited[n] = 1;
	      stack[end] = n;
	      end++;
	    }
	  }
	  
	}
      }
    }
  }
}

