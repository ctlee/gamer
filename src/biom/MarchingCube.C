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
 * File:     MarchingCube.C    < ... >
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Calculate the isosurface of a volume
 *           
 * Source:   Part of this file was adapted from an implementation by:           
 *           Copyright (c) 1998 Michael A. Chupa.
 *           Permission is expressly granted to use this code in any 
 *           non-commercial work, provided that this notice is preserved. 
 *           However, be aware that the Marching Cubes algorithm is patented.   
 * ***************************************************************************
 */


#include "MarchingCube.h"
 
SurfaceMesh* SurfaceMesh_marchingCube(int xdim, int ydim, int zdim, float* dataset, 
				      float isovalue, float* intensity, 
				      float intensity_isovalue, SPNT **holelist);

float get_intensity_ratio(float den1, float den2, float isovalue)
{
  float ratio;

  // If we are on the padded border
  if (den1 < 0 || den2 < 0)
    return 0.5;

  if (den1 != den2)
  {
    ratio = (isovalue-den1)/(den2-den1);
    if (ratio < 0.1)
      return 0.1;
    if (ratio > 0.9)
      return 0.9;
    return ratio;
  }
  return 0.5;
}

SurfaceMesh* SurfaceMesh_marchingCube(int xdim, int ydim, int zdim, 
				      float* dataset, float isovalue, 
				      SPNT **holelist)
{
  return SurfaceMesh_marchingCube(xdim, ydim, zdim, dataset, isovalue, 
				  NULL, 0.0, holelist);
}

SurfaceMesh* SurfaceMesh_marchingCube(int xdim, int ydim, int zdim, float* dataset, 
				      float isovalue, float* intensity, 
				      float intensity_isovalue, SPNT **holelist)
{
  int tempt_x, tempt_y, tempt_z;
  int i,j,k;
  int m,n,l;
  int number, ii, stack_size;
  float den1,den2,ratio;
  int v_num, t_num;
  FLTVECT *vertex;
  INT3VECT *triangle;
  int cellVerts[12]; 
  unsigned char cellIndex;      
  INT3VECT *mc_edge;
  unsigned char *mc_sign;
  SurfaceMesh *surfmesh;
  SPNT *hole_end, *hole_start;
  FLTVECT interior_seed;
  float max_density = -9999.0;
  bool use_intensity;

  // Check if we are using intensity values (neede for Lattice Data)
  use_intensity = intensity != NULL;

  // Initialize memory
  vertex = (FLTVECT*)malloc(sizeof(FLTVECT)*xdim*ydim*zdim);
  triangle = (INT3VECT*)malloc(sizeof(INT3VECT)*xdim*ydim*zdim);
  v_num = 0;
  t_num = 0;
  
  mc_edge = (INT3VECT*)malloc(sizeof(INT3VECT)*xdim*ydim*zdim);
  mc_sign = (unsigned char*)malloc(sizeof(unsigned char)*xdim*ydim*zdim);

  // preprocessing: removing small holes
  for (k=0; k<zdim; k++)
    for (j=0; j<ydim; j++) 
      for (i=0; i<xdim; i++) {
	if (dataset[IndexVect(i,j,k)] > max_density)
	  max_density = dataset[IndexVect(i,j,k)];
	mc_sign[IndexVect(i,j,k)] = 0;
      }

  triangle[0].a = 0;
  triangle[0].b = 0;
  triangle[0].c = 0;
  mc_sign[0] = 1;
  stack_size = 1;
  while (stack_size > 0) {
    stack_size--;
    tempt_x = triangle[stack_size].a;
    tempt_y = triangle[stack_size].b;
    tempt_z = triangle[stack_size].c;
    for (k=max(tempt_z-1,0); k<=min(tempt_z+1,zdim-1); k++) 
      for (j=max(tempt_y-1,0); j<=min(tempt_y+1,ydim-1); j++) 
	for (i=max(tempt_x-1,0); i<=min(tempt_x+1,xdim-1); i++) {
	  if (dataset[IndexVect(i,j,k)] < isovalue &&
	      mc_sign[IndexVect(i,j,k)] == 0) {
	    mc_sign[IndexVect(i,j,k)] = 1;
	    triangle[stack_size].a = i;
	    triangle[stack_size].b = j;
	    triangle[stack_size].c = k;
	    stack_size++;
	  }
	}
  }   
  hole_start = NULL;
  hole_end = NULL;
  den1 = 0;
  for (l=0; l<zdim; l++)
    for (n=0; n<ydim; n++) 
      for (m=0; m<xdim; m++) {
	if (dataset[IndexVect(m,n,l)] > den1) {
	  den1 = dataset[IndexVect(m,n,l)];
	  interior_seed.x = (float)m;
	  interior_seed.y = (float)n;
	  interior_seed.z = (float)l;
	}

	if (dataset[IndexVect(m,n,l)] < isovalue &&
	    mc_sign[IndexVect(m,n,l)] == 0) {
	  number = 1;
	  triangle[0].a = m;
	  triangle[0].b = n;
	  triangle[0].c = l;
	  mc_sign[0] = 1;
	  stack_size = 1;
	  while (stack_size > 0) {
	    stack_size--;
	    tempt_x = triangle[stack_size].a;
	    tempt_y = triangle[stack_size].b;
	    tempt_z = triangle[stack_size].c;
	    for (k=max(tempt_z-1,0); k<=min(tempt_z+1,zdim-1); k++) 
	      for (j=max(tempt_y-1,0); j<=min(tempt_y+1,ydim-1); j++) 
		for (i=max(tempt_x-1,0); i<=min(tempt_x+1,xdim-1); i++) {
		  if (dataset[IndexVect(i,j,k)] < isovalue &&
		      mc_sign[IndexVect(i,j,k)] == 0) {
		    mc_sign[IndexVect(i,j,k)] = 1;
		    triangle[stack_size].a = i;
		    triangle[stack_size].b = j;
		    triangle[stack_size].c = k;
		    stack_size++;
		    number++;
		  }
		}
	  }
	  printf("hole size: %d \n",number);
	  if (number < MIN_VOLUME) {
	    triangle[0].a = m;
	    triangle[0].b = n;
	    triangle[0].c = l;
	    dataset[IndexVect(m,n,l)] = max_density;
	    stack_size = 1;
	    while (stack_size > 0) {
	      stack_size--;
	      tempt_x = triangle[stack_size].a;
	      tempt_y = triangle[stack_size].b;
	      tempt_z = triangle[stack_size].c;
	      for (k=max(tempt_z-1,0); k<=min(tempt_z+1,zdim-1); k++) 
		for (j=max(tempt_y-1,0); j<=min(tempt_y+1,ydim-1); j++) 
		  for (i=max(tempt_x-1,0); i<=min(tempt_x+1,xdim-1); i++) {
		    if (dataset[IndexVect(i,j,k)] < isovalue) {
		      dataset[IndexVect(i,j,k)] = max_density;
		      triangle[stack_size].a = i;
		      triangle[stack_size].b = j;
		      triangle[stack_size].c = k;
		      stack_size++;
		    }
		  }
	    }
	  }
	  else {
	    /* could be improved here ?????*/ 
	    if (hole_start == NULL) {
	      hole_end = (SPNT*)malloc(sizeof(SPNT));
	      hole_start = hole_end;
	    }
	    else {
	      hole_end->next = (SPNT*)malloc(sizeof(SPNT));
	      hole_end = hole_end->next;
	    }
	    hole_end->x = (float)m;
	    hole_end->y = (float)n;
	    hole_end->z = (float)l;
	  }
	}
      }
  if (hole_end != NULL)
    hole_end->next = NULL;
  *holelist = hole_start;
  
  for (k=0; k<zdim; k++)
    for (j=0; j<ydim; j++) 
      for (i=0; i<xdim; i++) {
	if (dataset[IndexVect(i,j,k)] > isovalue-0.0001 &&
	    dataset[IndexVect(i,j,k)] < isovalue+0.0001)
	  dataset[IndexVect(i,j,k)] = isovalue+0.0001;
	
	mc_edge[IndexVect(i,j,k)].a = -1;
	mc_edge[IndexVect(i,j,k)].b = -1;
	mc_edge[IndexVect(i,j,k)].c = -1;
	if (dataset[IndexVect(i,j,k)] >= isovalue)
	  mc_sign[IndexVect(i,j,k)] = 1;
	else
	  mc_sign[IndexVect(i,j,k)] = 255;
      }	
  
  for (tempt_z=0; tempt_z<zdim-1; tempt_z++)
    for (tempt_y=0; tempt_y<ydim-1; tempt_y++) 
      for (tempt_x=0; tempt_x<xdim-1; tempt_x++) {
	
	for (ii = 0; ii < 12; ii++)
	  cellVerts[ii] = -1;
	
	cellIndex = 0;
	if (mc_sign[IndexVect(tempt_x,tempt_y,tempt_z)] == 255) cellIndex |= 1;
	if (mc_sign[IndexVect(tempt_x,tempt_y+1,tempt_z)] == 255) cellIndex |= 2;
	if (mc_sign[IndexVect(tempt_x+1,tempt_y+1,tempt_z)] == 255) cellIndex |= 4;
	if (mc_sign[IndexVect(tempt_x+1,tempt_y,tempt_z)] == 255) cellIndex |= 8;
	if (mc_sign[IndexVect(tempt_x,tempt_y,tempt_z+1)] == 255) cellIndex |= 16;
	if (mc_sign[IndexVect(tempt_x,tempt_y+1,tempt_z+1)] == 255) cellIndex |= 32;
	if (mc_sign[IndexVect(tempt_x+1,tempt_y+1,tempt_z+1)] == 255) cellIndex |= 64;
	if (mc_sign[IndexVect(tempt_x+1,tempt_y,tempt_z+1)] == 255) cellIndex |= 128;     
	
	if (edgeTable[cellIndex] & 1) {
	  if (mc_edge[IndexVect(tempt_x, tempt_y, tempt_z)].b == -1) {  
	    
	    if(use_intensity)
	    {
	      den1 = intensity[IndexVect(tempt_x, tempt_y, tempt_z)];
	      den2 = intensity[IndexVect(tempt_x, tempt_y+1, tempt_z)];

	      ratio = get_intensity_ratio(den1, den2, intensity_isovalue);
  	    }
	    else
	    {
	      den1 = dataset[IndexVect(tempt_x, tempt_y, tempt_z)];
	      den2 = dataset[IndexVect(tempt_x, tempt_y+1, tempt_z)];
	      if (den1 != den2)
		ratio = (isovalue-den1)/(den2-den1);
	      else 
		ratio = 0;
	    }
	    vertex[v_num].x = (float)tempt_x;
	    vertex[v_num].y = (float)tempt_y+ratio;
	    vertex[v_num].z = (float)tempt_z;
	    cellVerts[0] = v_num;
	    mc_edge[IndexVect(tempt_x,tempt_y,tempt_z)].b = v_num;
	    v_num++;
	  }
	  else 
	    cellVerts[0] = mc_edge[IndexVect(tempt_x, tempt_y, tempt_z)].b;
	}

	if (edgeTable[cellIndex] & 2) {
	  if (mc_edge[IndexVect(tempt_x, tempt_y+1, tempt_z)].a == -1) {  
	    if(use_intensity)
	    {
	      den1 = intensity[IndexVect(tempt_x, tempt_y+1, tempt_z)];
	      den2 = intensity[IndexVect(tempt_x+1, tempt_y+1, tempt_z)];

	      ratio = get_intensity_ratio(den1, den2, intensity_isovalue);
  	    }
	    else
	    {
	      den1 = dataset[IndexVect(tempt_x, tempt_y+1, tempt_z)];
	      den2 = dataset[IndexVect(tempt_x+1, tempt_y+1, tempt_z)];
	      if (den1 != den2)
		ratio = (isovalue-den1)/(den2-den1);
	      else 
		ratio = 0;
	    }
	    vertex[v_num].x = (float)tempt_x+ratio;
	    vertex[v_num].y = (float)tempt_y+1;
	    vertex[v_num].z = (float)tempt_z;
	    cellVerts[1]  = v_num;
	    mc_edge[IndexVect(tempt_x, tempt_y+1, tempt_z)].a = v_num;
	    v_num++;
	  }
	  else 
	    cellVerts[1] = mc_edge[IndexVect(tempt_x, tempt_y+1, tempt_z)].a;
	}
	if (edgeTable[cellIndex] & 4) {
	  if (mc_edge[IndexVect(tempt_x+1,tempt_y,tempt_z)].b == -1) {  
	    if(use_intensity)
	    {
	      den1 = intensity[IndexVect(tempt_x+1, tempt_y, tempt_z)];
	      den2 = intensity[IndexVect(tempt_x+1, tempt_y+1, tempt_z)];

	      ratio = get_intensity_ratio(den1, den2, intensity_isovalue);
  	    }
	    else
	    {
	      den1 = dataset[IndexVect(tempt_x+1, tempt_y, tempt_z)];
	      den2 = dataset[IndexVect(tempt_x+1, tempt_y+1, tempt_z)];
	      if (den1 != den2)
		ratio = (isovalue-den1)/(den2-den1);
	      else 
		ratio = 0;
	    }
	    vertex[v_num].x = (float)tempt_x+1;
	    vertex[v_num].y = (float)tempt_y+ratio;
	    vertex[v_num].z = (float)tempt_z;
	    cellVerts[2] = v_num;
	    mc_edge[IndexVect(tempt_x+1, tempt_y, tempt_z)].b = v_num;
	    v_num++;
	  }
	  else 
	    cellVerts[2] = mc_edge[IndexVect(tempt_x+1, tempt_y, tempt_z)].b;
	}
	if (edgeTable[cellIndex] & 8) {
	  if (mc_edge[IndexVect(tempt_x, tempt_y, tempt_z)].a == -1) {  
	    if(use_intensity)
	    {
	      den1 = intensity[IndexVect(tempt_x, tempt_y, tempt_z)];
	      den2 = intensity[IndexVect(tempt_x+1, tempt_y, tempt_z)];

	      ratio = get_intensity_ratio(den1, den2, intensity_isovalue);
  	    }
	    else
	    {
	      den1 = dataset[IndexVect(tempt_x, tempt_y, tempt_z)];
	      den2 = dataset[IndexVect(tempt_x+1, tempt_y, tempt_z)];
	      if (den1 != den2)
		ratio = (isovalue-den1)/(den2-den1);
	      else 
		ratio = 0;
	    }
	    vertex[v_num].x = (float)tempt_x+ratio;
	    vertex[v_num].y = (float)tempt_y;
	    vertex[v_num].z = (float)tempt_z;
	    cellVerts[3]  = v_num;
	    mc_edge[IndexVect(tempt_x, tempt_y, tempt_z)].a = v_num;
	    v_num++;
	  }
	  else 
	    cellVerts[3] = mc_edge[IndexVect(tempt_x, tempt_y, tempt_z)].a;
	}
	if (edgeTable[cellIndex] & 16) {
	  if (mc_edge[IndexVect(tempt_x,tempt_y,tempt_z+1)].b == -1) {  
	    if(use_intensity)
	    {
	      den1 = intensity[IndexVect(tempt_x, tempt_y, tempt_z+1)];
	      den2 = intensity[IndexVect(tempt_x, tempt_y+1, tempt_z+1)];

	      ratio = get_intensity_ratio(den1, den2, intensity_isovalue);
  	    }
	    else
	    {
	      den1 = dataset[IndexVect(tempt_x, tempt_y, tempt_z+1)];
	      den2 = dataset[IndexVect(tempt_x, tempt_y+1, tempt_z+1)];
	      if (den1 != den2)
		ratio = (isovalue-den1)/(den2-den1);
	      else 
		ratio = 0;
	    }
	    vertex[v_num].x = (float)tempt_x;
	    vertex[v_num].y = (float)tempt_y+ratio;
	    vertex[v_num].z = (float)tempt_z+1;
	    cellVerts[4]  = v_num;
	    mc_edge[IndexVect(tempt_x, tempt_y, tempt_z+1)].b = v_num;
	    v_num++;
	  }
	  else 
	    cellVerts[4] = mc_edge[IndexVect(tempt_x, tempt_y, tempt_z+1)].b;
	}
	if (edgeTable[cellIndex] & 32) {
	  if (mc_edge[IndexVect(tempt_x,tempt_y+1,tempt_z+1)].a == -1) {  
	    if(use_intensity)
	    {
	      den1 = intensity[IndexVect(tempt_x, tempt_y+1, tempt_z+1)];
	      den2 = intensity[IndexVect(tempt_x+1, tempt_y+1, tempt_z+1)];

	      ratio = get_intensity_ratio(den1, den2, intensity_isovalue);
  	    }
	    else
	    {
	      den1 = dataset[IndexVect(tempt_x, tempt_y+1, tempt_z+1)];
	      den2 = dataset[IndexVect(tempt_x+1, tempt_y+1, tempt_z+1)];
	      if (den1 != den2)
		ratio = (isovalue-den1)/(den2-den1);
	      else 
		ratio = 0;
	    }
	    vertex[v_num].x = (float)tempt_x+ratio;
	    vertex[v_num].y = (float)tempt_y+1;
	    vertex[v_num].z = (float)tempt_z+1;
	    cellVerts[5]  = v_num;
	    mc_edge[IndexVect(tempt_x, tempt_y+1, tempt_z+1)].a = v_num;
	    v_num++;
	  }
	  else 
	    cellVerts[5] = mc_edge[IndexVect(tempt_x, tempt_y+1, tempt_z+1)].a;
	}
	if (edgeTable[cellIndex] & 64) {
	  if (mc_edge[IndexVect(tempt_x+1, tempt_y, tempt_z+1)].b == -1) {  
	    if(use_intensity)
	    {
	      den1 = intensity[IndexVect(tempt_x+1, tempt_y, tempt_z+1)];
	      den2 = intensity[IndexVect(tempt_x+1, tempt_y+1, tempt_z+1)];

	      ratio = get_intensity_ratio(den1, den2, intensity_isovalue);
  	    }
	    else
	    {
	      den1 = dataset[IndexVect(tempt_x+1, tempt_y, tempt_z+1)];
	      den2 = dataset[IndexVect(tempt_x+1, tempt_y+1, tempt_z+1)];
	      if (den1 != den2)
		ratio = (isovalue-den1)/(den2-den1);
	      else 
		ratio = 0;
	    }
	    vertex[v_num].x = (float)tempt_x+1;
	    vertex[v_num].y = (float)tempt_y+ratio;
	    vertex[v_num].z = (float)tempt_z+1;
	    cellVerts[6]  = v_num;
	    mc_edge[IndexVect(tempt_x+1, tempt_y, tempt_z+1)].b = v_num;
	    v_num++;
	  }
	  else 
	    cellVerts[6] = mc_edge[IndexVect(tempt_x+1, tempt_y, tempt_z+1)].b;
	}
	if (edgeTable[cellIndex] & 128) {
	  if (mc_edge[IndexVect(tempt_x, tempt_y, tempt_z+1)].a == -1) {  
	    if(use_intensity)
	    {
	      den1 = intensity[IndexVect(tempt_x, tempt_y, tempt_z+1)];
	      den2 = intensity[IndexVect(tempt_x+1, tempt_y, tempt_z+1)];

	      ratio = get_intensity_ratio(den1, den2, intensity_isovalue);
  	    }
	    else
	    {
	      den1 = dataset[IndexVect(tempt_x, tempt_y, tempt_z+1)];
	      den2 = dataset[IndexVect(tempt_x+1, tempt_y, tempt_z+1)];
	      if (den1 != den2)
		ratio = (isovalue-den1)/(den2-den1);
	      else 
		ratio = 0;
	    }
	    vertex[v_num].x = (float)tempt_x+ratio;
	    vertex[v_num].y = (float)tempt_y;
	    vertex[v_num].z = (float)tempt_z+1;
	    cellVerts[7]  = v_num;
	    mc_edge[IndexVect(tempt_x, tempt_y, tempt_z+1)].a = v_num;
	    v_num++;
	  }
	  else 
	    cellVerts[7] = mc_edge[IndexVect(tempt_x, tempt_y, tempt_z+1)].a;
	}
	if (edgeTable[cellIndex] & 256) {
	  if (mc_edge[IndexVect(tempt_x, tempt_y, tempt_z)].c == -1) {  
	    if(use_intensity)
	    {
	      den1 = intensity[IndexVect(tempt_x, tempt_y, tempt_z)];
	      den2 = intensity[IndexVect(tempt_x, tempt_y, tempt_z+1)];

	      ratio = get_intensity_ratio(den1, den2, intensity_isovalue);
  	    }
	    else
	    {
	      den1 = dataset[IndexVect(tempt_x, tempt_y, tempt_z)];
	      den2 = dataset[IndexVect(tempt_x, tempt_y, tempt_z+1)];
	      if (den1 != den2)
		ratio = (isovalue-den1)/(den2-den1);
	      else 
		ratio = 0;
	    }
	    vertex[v_num].x = (float)tempt_x;
	    vertex[v_num].y = (float)tempt_y;
	    vertex[v_num].z = (float)tempt_z+ratio;
	    cellVerts[8]  = v_num;
	    mc_edge[IndexVect(tempt_x, tempt_y, tempt_z)].c = v_num;
	    v_num++;
	  }
	  else 
	    cellVerts[8] = mc_edge[IndexVect(tempt_x,tempt_y,tempt_z)].c;
	}
	if (edgeTable[cellIndex] & 512) {
	  if (mc_edge[IndexVect(tempt_x, tempt_y+1, tempt_z)].c == -1) {  
	    if(use_intensity)
	    {
	      den1 = intensity[IndexVect(tempt_x, tempt_y+1, tempt_z)];
	      den2 = intensity[IndexVect(tempt_x, tempt_y+1, tempt_z+1)];

	      ratio = get_intensity_ratio(den1, den2, intensity_isovalue);
  	    }
	    else
	    {
	      den1 = dataset[IndexVect(tempt_x, tempt_y+1, tempt_z)];
	      den2 = dataset[IndexVect(tempt_x, tempt_y+1, tempt_z+1)];
	      if (den1 != den2)
		ratio = (isovalue-den1)/(den2-den1);
	      else 
		ratio = 0;
	    }
	    vertex[v_num].x = (float)tempt_x;
	    vertex[v_num].y = (float)tempt_y+1;
	    vertex[v_num].z = (float)tempt_z+ratio;
	    cellVerts[9]  = v_num;
	    mc_edge[IndexVect(tempt_x, tempt_y+1, tempt_z)].c = v_num;
	    v_num++;
	  }
	  else 
	    cellVerts[9] = mc_edge[IndexVect(tempt_x, tempt_y+1, tempt_z)].c;
	}
	if (edgeTable[cellIndex] & 1024) {
	  if (mc_edge[IndexVect(tempt_x+1, tempt_y+1, tempt_z)].c == -1) {  
	    if(use_intensity)
	    {
	      den1 = intensity[IndexVect(tempt_x+1, tempt_y+1, tempt_z)];
	      den2 = intensity[IndexVect(tempt_x+1, tempt_y+1, tempt_z+1)];

	      ratio = get_intensity_ratio(den1, den2, intensity_isovalue);
  	    }
	    else
	    {
	      den1 = dataset[IndexVect(tempt_x+1, tempt_y+1, tempt_z)];
	      den2 = dataset[IndexVect(tempt_x+1, tempt_y+1, tempt_z+1)];
	      if (den1 != den2)
		ratio = (isovalue-den1)/(den2-den1);
	      else 
		ratio = 0;
	    }
	    vertex[v_num].x = (float)tempt_x+1;
	    vertex[v_num].y = (float)tempt_y+1;
	    vertex[v_num].z = (float)tempt_z+ratio;
	    cellVerts[10]  = v_num;
	    mc_edge[IndexVect(tempt_x+1, tempt_y+1, tempt_z)].c = v_num;
	    v_num++;
	  }
	  else 
	    cellVerts[10] = mc_edge[IndexVect(tempt_x+1, tempt_y+1, tempt_z)].c;
	}
	if (edgeTable[cellIndex] & 2048) {
	  if (mc_edge[IndexVect(tempt_x+1,tempt_y,tempt_z)].c == -1) {  
	    if(use_intensity)
	    {
	      den1 = intensity[IndexVect(tempt_x+1, tempt_y, tempt_z)];
	      den2 = intensity[IndexVect(tempt_x+1, tempt_y, tempt_z+1)];

	      ratio = get_intensity_ratio(den1, den2, intensity_isovalue);
  	    }
	    else
	    {
	      den1 = dataset[IndexVect(tempt_x+1, tempt_y, tempt_z)];
	      den2 = dataset[IndexVect(tempt_x+1, tempt_y, tempt_z+1)];
	      if (den1 != den2)
		ratio = (isovalue-den1)/(den2-den1);
	      else 
		ratio = 0;
	    }
	    vertex[v_num].x = (float)tempt_x+1;
	    vertex[v_num].y = (float)tempt_y;
	    vertex[v_num].z = (float)tempt_z+ratio;
	    cellVerts[11]  = v_num;
	    mc_edge[IndexVect(tempt_x+1, tempt_y, tempt_z)].c = v_num;
	    v_num++;
	  }
	  else 
	    cellVerts[11] = mc_edge[IndexVect(tempt_x+1, tempt_y, tempt_z)].c;
	}
	
	ii = 0;
	while (triTable[cellIndex][ii] != -1) {
	  triangle[t_num].a = cellVerts[triTable[cellIndex][ii++]];
	  triangle[t_num].b = cellVerts[triTable[cellIndex][ii++]];
	  triangle[t_num].c = cellVerts[triTable[cellIndex][ii++]];
	  t_num++;
	}
      }

  // Allocate memory
  surfmesh = SurfaceMesh_ctor(v_num, t_num);
  
  for (n = 0; n < surfmesh->num_vertices; n++) {
    surfmesh->vertex[n].x = vertex[n].x;
    surfmesh->vertex[n].y = vertex[n].y;
    surfmesh->vertex[n].z = vertex[n].z;
  }
  for (n = 0; n < surfmesh->num_faces; n++) {
    surfmesh->face[n].a = triangle[n].a;
    surfmesh->face[n].b = triangle[n].b;
    surfmesh->face[n].c = triangle[n].c;
  }
  
  free(vertex);
  free(triangle);
  free(mc_edge);
  free(mc_sign);

  // Return created SurfaceMesh
  return surfmesh;
  
}

