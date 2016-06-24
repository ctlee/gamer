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
 * File:     GenerateSphereMesh.C    < ... >
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Generate sphere mesh used as a bounding mesh for molecules
 *
 * Notes:    Starting from an octahedron (6 vertices, 8 faces, and 12 edges)
 *           Consider level n (original level = 1):
 *           We can see that f_n = 8 * 4^(n-1),
 *           By recursive scheme: v_n = 2v_{n-1} + 8*4^{n-2} - 2, we can get:
 *           Number of vertices = 12*2^{n-2}+8*2^{2n-3}-8*2^{n-2}-2^n+2.
 *
 *           The following are the calculated v_n and f_n for n = 1, ... 9:
 *           level 1: v_n = 6,  f_n = 8 
 *           level 2: v_n = 18,  f_n = 32 
 *           level 3: v_n = 66,  f_n = 128 
 *           level 4: v_n = 258,  f_n = 512 
 *           level 5: v_n = 1026,  f_n = 2048 
 *           level 6: v_n = 4098,  f_n = 8192 
 *           level 7: v_n = 16386,  f_n = 32768 
 *           level 8: v_n = 65538,  f_n = 131072 
 *           level 9: v_n = 262146,  f_n = 524288 
 *           In most cases, n=4 should be a good choice. try not to use over 7
 * ****************************************************************************
 */


#include <gamer/biom.h>
#include "gamercf.h"

SurfaceMesh* SurfaceMesh_sphere(int level)
{
  int n,m,l,k;
  int a,b,c,e,f,g;
  float ax,ay,az;
  float bx,by,bz;
  float cx,cy,cz;
  float ex,ey,ez;
  float fx,fy,fz;
  float gx,gy,gz;
  SurfaceMesh *surfmesh;
  int n_v, n_f;
  float length,dist;

  // Allocate memory
  if (level == 1)
    surfmesh = SurfaceMesh_ctor(6, 8*(int)pow(4,(level-1)));
  else if (level > 1)
    surfmesh = SurfaceMesh_ctor(12*(int)pow(2,(level-2)) + 8*(int)pow(2,(2*level-3)) - 
				8*(int)pow(2,(level-2)) - (int)pow(2,level) + 2, 
				8*(int)pow(4,(level-1)));
  else {
    printf("  The level of subdivision must be greater than 0 (usually 4 is good) ....\n");
    printf("  Try not to use over 7, otherwise, you may have to wait for minutes/hours !!!! .... \n");
    exit(0);
  }

  printf("vertices: %d --- triangles: %d \n",surfmesh->num_vertices,surfmesh->num_faces);
    
  // initialization
  n_v = 6;
  n_f = 8;
  surfmesh->vertex[0].x = 0;
  surfmesh->vertex[0].y = 0;
  surfmesh->vertex[0].z = 1;
  surfmesh->vertex[1].x = 0;
  surfmesh->vertex[1].y = 0;
  surfmesh->vertex[1].z = -1;
  surfmesh->vertex[2].x = 0;
  surfmesh->vertex[2].y = 1;
  surfmesh->vertex[2].z = 0;
  surfmesh->vertex[3].x = 0;
  surfmesh->vertex[3].y = -1;
  surfmesh->vertex[3].z = 0;
  surfmesh->vertex[4].x = 1;
  surfmesh->vertex[4].y = 0;
  surfmesh->vertex[4].z = 0;
  surfmesh->vertex[5].x = -1;
  surfmesh->vertex[5].y = 0;
  surfmesh->vertex[5].z = 0;
  surfmesh->face[0].a = 0;
  surfmesh->face[0].b = 4;
  surfmesh->face[0].c = 2;
  surfmesh->face[1].a = 0;
  surfmesh->face[1].b = 2;
  surfmesh->face[1].c = 5;
  surfmesh->face[2].a = 0;
  surfmesh->face[2].b = 5;
  surfmesh->face[2].c = 3;
  surfmesh->face[3].a = 0;
  surfmesh->face[3].b = 3;
  surfmesh->face[3].c = 4;
  surfmesh->face[4].a = 1;
  surfmesh->face[4].b = 2;
  surfmesh->face[4].c = 4;
  surfmesh->face[5].a = 1;
  surfmesh->face[5].b = 5;
  surfmesh->face[5].c = 2;
  surfmesh->face[6].a = 1;
  surfmesh->face[6].b = 3;
  surfmesh->face[6].c = 5;
  surfmesh->face[7].a = 1;
  surfmesh->face[7].b = 4;
  surfmesh->face[7].c = 3;
  
  // subdivision
  m = 0;
  while (m < level-1) {
    surfmesh->num_vertices = n_v;
    surfmesh->num_faces = n_f;
    for (n = 0; n < surfmesh->num_faces; n++) {
      a = surfmesh->face[n].a;
      b = surfmesh->face[n].b;
      c = surfmesh->face[n].c;
      ax = surfmesh->vertex[a].x;
      ay = surfmesh->vertex[a].y;
      az = surfmesh->vertex[a].z;
      bx = surfmesh->vertex[b].x;
      by = surfmesh->vertex[b].y;
      bz = surfmesh->vertex[b].z;
      cx = surfmesh->vertex[c].x;
      cy = surfmesh->vertex[c].y;
      cz = surfmesh->vertex[c].z;
      
      ex = 0.5*(ax+bx);
      ey = 0.5*(ay+by);
      ez = 0.5*(az+bz);
      length = sqrt(ex*ex+ey*ey+ez*ez);
      ex /= length;
      ey /= length;
      ez /= length;
      k = -1;
      for (l = 0; l < n_v; l++) {
	dist = sqrt((ex-surfmesh->vertex[l].x)*(ex-surfmesh->vertex[l].x)+
		    (ey-surfmesh->vertex[l].y)*(ey-surfmesh->vertex[l].y)+
		    (ez-surfmesh->vertex[l].z)*(ez-surfmesh->vertex[l].z));
	if (dist < 0.0001) {
	  k = l;
	  break;
	}
      }//check if the vertex has been created or not
      if (k == -1) {
	e = n_v;
	n_v++;
      }
      else
	e = k;
      surfmesh->vertex[e].x = ex;
      surfmesh->vertex[e].y = ey;
      surfmesh->vertex[e].z = ez;
      
      fx = 0.5*(cx+bx);
      fy = 0.5*(cy+by);
      fz = 0.5*(cz+bz);
      length = sqrt(fx*fx+fy*fy+fz*fz);
      fx /= length;
      fy /= length;
      fz /= length;
      k = -1;
      for (l = 0; l < n_v; l++) {
	dist = sqrt((fx-surfmesh->vertex[l].x)*(fx-surfmesh->vertex[l].x)+
		    (fy-surfmesh->vertex[l].y)*(fy-surfmesh->vertex[l].y)+
		    (fz-surfmesh->vertex[l].z)*(fz-surfmesh->vertex[l].z));
	if (dist < 0.0001) {
	  k = l;
	  break;
	}
      }//check if the vertex has been created or not
      if (k == -1) {
	f = n_v;
	n_v++;
      }
      else
	f = k;
      surfmesh->vertex[f].x = fx;
      surfmesh->vertex[f].y = fy;
      surfmesh->vertex[f].z = fz;
      
      gx = 0.5*(cx+ax);
      gy = 0.5*(cy+ay);
      gz = 0.5*(cz+az);
      length = sqrt(gx*gx+gy*gy+gz*gz);
      gx /= length;
      gy /= length;
      gz /= length;
      k = -1;
      for (l = 0; l < n_v; l++) {
	dist = sqrt((gx-surfmesh->vertex[l].x)*(gx-surfmesh->vertex[l].x)+
		    (gy-surfmesh->vertex[l].y)*(gy-surfmesh->vertex[l].y)+
		    (gz-surfmesh->vertex[l].z)*(gz-surfmesh->vertex[l].z));
	if (dist < 0.0001) {
	  k = l;
	  break;
	}
      }//check if the vertex has been created or not
      if (k == -1) {
	g = n_v;
	n_v++;
      }
      else
	g = k;
      surfmesh->vertex[g].x = gx;
      surfmesh->vertex[g].y = gy;
      surfmesh->vertex[g].z = gz;
      
      // create new triangles
      surfmesh->face[n].a = e;
      surfmesh->face[n].b = f;
      surfmesh->face[n].c = g;
      surfmesh->face[n_f].a = a;
      surfmesh->face[n_f].b = e;
      surfmesh->face[n_f].c = g;
      n_f++;
      surfmesh->face[n_f].a = e;
      surfmesh->face[n_f].b = b;
      surfmesh->face[n_f].c = f;
      n_f++;
      surfmesh->face[n_f].a = g;
      surfmesh->face[n_f].b = f;
      surfmesh->face[n_f].c = c;
      n_f++;
    }
    m++;
  }
  surfmesh->num_vertices = n_v;
  surfmesh->num_faces = n_f;

  return surfmesh;
}
