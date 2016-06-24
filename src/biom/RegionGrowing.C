/*
 * ***************************************************************************
 * BIMoS = < Biomedical Image-based Modeling and Simulation >
 * Copyright (C) 2009-2010 -- Zeyun Yu (yuz@uwm.edu)
 * Dept. Computer Science, The University of Wisconsin-Milwaukee
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
 * File:     RegionGrowing.C    < ... >
 *
 * Author:   Zeyun Yu 
 *
 * Purpose:  Find the SES using the fast marching method
 * ***************************************************************************
 */


#include <gamer/biom.h>
#include "gamercf.h"


#define MaxDist    29999

int xdim1,ydim1,zdim1;
MinHeapS* min_heap;
SEEDS *AllSeeds;
int *heap_pointer;
ATOM *atom_list;
float threshold;
unsigned short min_x,min_y,min_z;
int min_seed;
float min_dist;


void GetMinimum(void);
void InsertHeap(int, int, int, float); 
void UpdateHeap(int, int, int,float);
void Marching(void);
FLTVECT FindSeed(float, float, float, int);



void ExtractSES(MinHeapS *mheap, SEEDS *all_seeds, int *heappointer, int xd, int yd, int zd,
		int *atom_index, int atom_num, ATOM *atomlist, float thresh)
{
  int i,j,k;
  int m,n,l, num, c;
  int index,index1;
  float dist;
  FLTVECT seed;
  char visited;


  xdim1 = xd;
  ydim1 = yd;
  zdim1 = zd;
  atom_list = atomlist;
  min_heap = mheap;
  AllSeeds = all_seeds;
  heap_pointer = heappointer;
  threshold = thresh;

  /* Initialize */ 
  index = 0;
  min_heap->size = 0;
  for (k=0; k<zdim1; k++)
    for (j=0; j<ydim1; j++)
      for (i=0; i<xdim1; i++) {
	if (atom_index[IndexVect1(i,j,k)] < 0) { 
	  for (num = 0; num < MaxAtom; num++)
	    AllSeeds[index].atom[num] = -1;
	  num = 0;
	  for (l=k-1; l<=k+1; l++) 
	    for (n=j-1; n<=j+1; n++) 
	      for (m=i-1; m<=i+1; m++) {
		if (m==i||n==j||l==k) {
		  index1 = atom_index[IndexVect1(m,n,l)];
		  if (index1 < 0) {
		    index1 = -index1-1;
		    visited = 0;
		    for (c=0; c<num; c++) {
		      if (index1 == AllSeeds[index].atom[c])
			visited = 1;
		    }
		    if (visited == 0) {
		      AllSeeds[index].atom[num] = index1;		   
		      num++;
		      if (num == MaxAtom)
			num--;
		    }
		  }
		}
	      }
	  	  
	  seed = FindSeed(i,j,k,index);
	  AllSeeds[index].seedx = seed.x;
	  AllSeeds[index].seedy = seed.y;
	  AllSeeds[index].seedz = seed.z;
	  dist = (seed.x-i)*(seed.x-i) + (seed.y-j)*(seed.y-j) + (seed.z-k)*(seed.z-k);
	  min_seed = index;
	  InsertHeap(i,j,k, dist);
	
	  index++;
	}
	else if (atom_index[IndexVect1(i,j,k)] > 0) {
	  heap_pointer[IndexVect1(i,j,k)] = MaxVal;
	}
	else {
	  heap_pointer[IndexVect1(i,j,k)] = -11;
        }
      }


  /* Fast Marching Method */
  while (1){
    GetMinimum();
    if (min_dist >= MaxDist-0.001)
      break;

    Marching();    
  }

}




void GetMinimum(void)
{
  int pointer, left, right;
  float dist;

  min_x = min_heap->x[0];
  min_y = min_heap->y[0];
  min_z = min_heap->z[0];
  min_seed = min_heap->seed[0];
  min_dist = min_heap->dist[0];


  if (min_dist == MaxDist)
    return;

  heap_pointer[IndexVect1(min_heap->x[0],min_heap->y[0],min_heap->z[0])] = -3;

  min_heap->size--;
  dist=min_heap->dist[min_heap->size];
  
  pointer=1;
  while (pointer <= min_heap->size/2) {
    left=2*pointer;
    right=2*pointer+1;
    if ((min_heap->dist[left-1] <= min_heap->dist[right-1]) && (min_heap->dist[left-1] < dist)) {
      min_heap->x[pointer-1]=min_heap->x[left-1];
      min_heap->y[pointer-1]=min_heap->y[left-1];
      min_heap->z[pointer-1]=min_heap->z[left-1];
      min_heap->seed[pointer-1]=min_heap->seed[left-1];
      min_heap->dist[pointer-1]=min_heap->dist[left-1];
      heap_pointer[IndexVect1(min_heap->x[pointer-1],min_heap->y[pointer-1],min_heap->z[pointer-1])] = pointer-1;
      pointer=left;
    }
    else if ((min_heap->dist[left-1] > min_heap->dist[right-1]) && (min_heap->dist[right-1] < dist)){
      min_heap->x[pointer-1]=min_heap->x[right-1];
      min_heap->y[pointer-1]=min_heap->y[right-1];
      min_heap->z[pointer-1]=min_heap->z[right-1];
      min_heap->seed[pointer-1]=min_heap->seed[right-1];
      min_heap->dist[pointer-1]=min_heap->dist[right-1];
      heap_pointer[IndexVect1(min_heap->x[pointer-1],min_heap->y[pointer-1],min_heap->z[pointer-1])] = pointer-1;
      pointer=right;
    }
    else break;
  }

  min_heap->x[pointer-1]=min_heap->x[min_heap->size];
  min_heap->y[pointer-1]=min_heap->y[min_heap->size];
  min_heap->z[pointer-1]=min_heap->z[min_heap->size];
  min_heap->seed[pointer-1]=min_heap->seed[min_heap->size];
  min_heap->dist[pointer-1]=dist;
  heap_pointer[IndexVect1(min_heap->x[min_heap->size],min_heap->y[min_heap->size],min_heap->z[min_heap->size])] = pointer-1;

}
    
void InsertHeap(int x, int y, int z, float dist)
{
  int pointer, parent;
  

  min_heap->size++;
  pointer=min_heap->size;

  while (pointer > 1) {
    if (pointer%2 == 0) {
      parent=pointer/2;
      if (dist < min_heap->dist[parent-1]) {
	min_heap->x[pointer-1]=min_heap->x[parent-1];
	min_heap->y[pointer-1]=min_heap->y[parent-1];
	min_heap->z[pointer-1]=min_heap->z[parent-1];
	min_heap->seed[pointer-1]=min_heap->seed[parent-1];
	min_heap->dist[pointer-1]=min_heap->dist[parent-1];
	heap_pointer[IndexVect1(min_heap->x[pointer-1],min_heap->y[pointer-1],min_heap->z[pointer-1])] = pointer-1;
	pointer=parent;
      }
      else break;
    }
    else if (pointer%2 == 1){
      parent=(pointer-1)/2;
      if (dist < min_heap->dist[parent-1]) {
	min_heap->x[pointer-1]=min_heap->x[parent-1];
	min_heap->y[pointer-1]=min_heap->y[parent-1];
	min_heap->z[pointer-1]=min_heap->z[parent-1];
	min_heap->seed[pointer-1]=min_heap->seed[parent-1];
	min_heap->dist[pointer-1]=min_heap->dist[parent-1];
	heap_pointer[IndexVect1(min_heap->x[pointer-1],min_heap->y[pointer-1],min_heap->z[pointer-1])] = pointer-1;
	pointer=parent;
      }
      else break;
    }
  }
  min_heap->x[pointer-1]=x;
  min_heap->y[pointer-1]=y;
  min_heap->z[pointer-1]=z;
  min_heap->seed[pointer-1]=min_seed;
  min_heap->dist[pointer-1]=dist;

  heap_pointer[IndexVect1(x,y,z)] = pointer-1;

}


void UpdateHeap(int x, int y, int z, float dist)
{
  int pointer, parent;
  int left, right;
  char up = 0;
  

  pointer=heap_pointer[IndexVect1(x,y,z)]+1;

  // checking the upper elements
  while (pointer > 1) {
    if (pointer%2 == 0) {
      parent=pointer/2;
      if (dist < min_heap->dist[parent-1]) {
	up = 1;
	min_heap->x[pointer-1]=min_heap->x[parent-1];
	min_heap->y[pointer-1]=min_heap->y[parent-1];
	min_heap->z[pointer-1]=min_heap->z[parent-1];
	min_heap->seed[pointer-1]=min_heap->seed[parent-1];
	min_heap->dist[pointer-1]=min_heap->dist[parent-1];
	heap_pointer[IndexVect1(min_heap->x[pointer-1],min_heap->y[pointer-1],min_heap->z[pointer-1])] = pointer-1; 
	pointer=parent;
      }
      else break;
    }
    else if (pointer%2 == 1){
      parent=(pointer-1)/2;
      if (dist < min_heap->dist[parent-1]) {
	up = 1;
	min_heap->x[pointer-1]=min_heap->x[parent-1];
	min_heap->y[pointer-1]=min_heap->y[parent-1];
	min_heap->z[pointer-1]=min_heap->z[parent-1];
	min_heap->seed[pointer-1]=min_heap->seed[parent-1];
	min_heap->dist[pointer-1]=min_heap->dist[parent-1];
	heap_pointer[IndexVect1(min_heap->x[pointer-1],min_heap->y[pointer-1],min_heap->z[pointer-1])] = pointer-1; 
	pointer=parent;
      }
      else break;
    }
  }

  if (up == 0) {
    // checking the lower elements
    while (pointer <= min_heap->size/2) {
      left=2*pointer;
      right=2*pointer+1;
      if ((min_heap->dist[left-1] <= min_heap->dist[right-1]) && (min_heap->dist[left-1] < dist)) {
	min_heap->x[pointer-1]=min_heap->x[left-1];
	min_heap->y[pointer-1]=min_heap->y[left-1];
	min_heap->z[pointer-1]=min_heap->z[left-1];
	min_heap->seed[pointer-1]=min_heap->seed[left-1];
	min_heap->dist[pointer-1]=min_heap->dist[left-1];
	heap_pointer[IndexVect1(min_heap->x[pointer-1],min_heap->y[pointer-1],min_heap->z[pointer-1])] = pointer-1;
	pointer=left;
      }
      else if ((min_heap->dist[left-1] > min_heap->dist[right-1]) && (min_heap->dist[right-1] < dist)){
	min_heap->x[pointer-1]=min_heap->x[right-1];
	min_heap->y[pointer-1]=min_heap->y[right-1];
	min_heap->z[pointer-1]=min_heap->z[right-1];
	min_heap->seed[pointer-1]=min_heap->seed[right-1];
	min_heap->dist[pointer-1]=min_heap->dist[right-1];
	heap_pointer[IndexVect1(min_heap->x[pointer-1],min_heap->y[pointer-1],min_heap->z[pointer-1])] = pointer-1;
	pointer=right;
      }
      else break;
    }
  }


  min_heap->x[pointer-1]=x;
  min_heap->y[pointer-1]=y;
  min_heap->z[pointer-1]=z;
  min_heap->seed[pointer-1]=min_seed;
  min_heap->dist[pointer-1]=dist;

  heap_pointer[IndexVect1(x,y,z)] = pointer-1;

}



void Marching(void)
{
  int tempt_x, tempt_y, tempt_z;
  float dt,dist;
  int neighbor,seed;
  char boundary;
  float min_seedx, min_seedy, min_seedz;
  float seedx, seedy, seedz;


  min_seedx = AllSeeds[min_seed].seedx;
  min_seedy = AllSeeds[min_seed].seedy;
  min_seedz = AllSeeds[min_seed].seedz;
  
  boundary = 0;
  

  /* ============================ */
  tempt_x=max(min_x-1,0);
  tempt_y=min_y;
  tempt_z=min_z;
  if (heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)] == MaxVal) {
    dist = (tempt_x-min_seedx)*(tempt_x-min_seedx)+(tempt_y-min_seedy)*(tempt_y-min_seedy)+(tempt_z-min_seedz)*(tempt_z-min_seedz);
    if (dist <= threshold)
      InsertHeap(tempt_x, tempt_y, tempt_z, dist);
    else 
      boundary = 1;
  }
  else if (heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)] > -1) {
    neighbor = heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)];
    dt = min_heap->dist[neighbor];
    if (dt < MaxDist) {
      dist = (tempt_x-min_seedx)*(tempt_x-min_seedx) + (tempt_y-min_seedy)*(tempt_y-min_seedy) + (tempt_z-min_seedz)*(tempt_z-min_seedz);
      if (dist < dt)
	UpdateHeap(tempt_x, tempt_y, tempt_z, dist);
    }
    else {
      dist = (tempt_x-min_seedx)*(tempt_x-min_seedx) + (tempt_y-min_seedy)*(tempt_y-min_seedy) + (tempt_z-min_seedz)*(tempt_z-min_seedz);
      seed = min_heap->seed[neighbor];
      seedx = AllSeeds[seed].seedx;
      seedy = AllSeeds[seed].seedy;
      seedz = AllSeeds[seed].seedz;
      if (dist < (tempt_x-seedx)*(tempt_x-seedx)+(tempt_y-seedy)*(tempt_y-seedy)+(tempt_z-seedz)*(tempt_z-seedz))
	UpdateHeap(tempt_x, tempt_y, tempt_z, MaxDist);
    }
  }
  //else if (heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)] > -11) {
      
  
  /* ============================ */
  tempt_x=min(min_x+1,xdim1-1);
  tempt_y=min_y;
  tempt_z=min_z;
  if (heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)] == MaxVal) {
    dist = (tempt_x-min_seedx)*(tempt_x-min_seedx)+(tempt_y-min_seedy)*(tempt_y-min_seedy)+(tempt_z-min_seedz)*(tempt_z-min_seedz);
    if (dist <= threshold)
      InsertHeap(tempt_x, tempt_y, tempt_z, dist);
    else 
      boundary = 1;
  }
  else if (heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)] > -1) {
    neighbor = heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)];
    dt = min_heap->dist[neighbor];
    if (dt < MaxDist) {
      dist = (tempt_x-min_seedx)*(tempt_x-min_seedx) + (tempt_y-min_seedy)*(tempt_y-min_seedy) + (tempt_z-min_seedz)*(tempt_z-min_seedz);
      if (dist < dt)
	UpdateHeap(tempt_x, tempt_y, tempt_z, dist);
    }
    else {
      dist = (tempt_x-min_seedx)*(tempt_x-min_seedx) + (tempt_y-min_seedy)*(tempt_y-min_seedy) + (tempt_z-min_seedz)*(tempt_z-min_seedz);
      seed = min_heap->seed[neighbor];
      seedx = AllSeeds[seed].seedx;
      seedy = AllSeeds[seed].seedy;
      seedz = AllSeeds[seed].seedz;
      if (dist < (tempt_x-seedx)*(tempt_x-seedx)+(tempt_y-seedy)*(tempt_y-seedy)+(tempt_z-seedz)*(tempt_z-seedz))
	UpdateHeap(tempt_x, tempt_y, tempt_z, MaxDist);
    }
  }
  

  /* ============================ */
  tempt_y=max(min_y-1,0);
  tempt_x=min_x;
  tempt_z=min_z;
  if (heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)] == MaxVal) {
    dist = (tempt_x-min_seedx)*(tempt_x-min_seedx)+(tempt_y-min_seedy)*(tempt_y-min_seedy)+(tempt_z-min_seedz)*(tempt_z-min_seedz);
    if (dist <= threshold)
      InsertHeap(tempt_x, tempt_y, tempt_z, dist);
    else 
      boundary = 1;
  }
  else if (heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)] > -1) {
    neighbor = heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)];
    dt = min_heap->dist[neighbor];
    if (dt < MaxDist) {
      dist = (tempt_x-min_seedx)*(tempt_x-min_seedx) + (tempt_y-min_seedy)*(tempt_y-min_seedy) + (tempt_z-min_seedz)*(tempt_z-min_seedz);
      if (dist < dt)
	UpdateHeap(tempt_x, tempt_y, tempt_z, dist);
    }
    else {
      dist = (tempt_x-min_seedx)*(tempt_x-min_seedx) + (tempt_y-min_seedy)*(tempt_y-min_seedy) + (tempt_z-min_seedz)*(tempt_z-min_seedz);
      seed = min_heap->seed[neighbor];
      seedx = AllSeeds[seed].seedx;
      seedy = AllSeeds[seed].seedy;
      seedz = AllSeeds[seed].seedz;
      if (dist < (tempt_x-seedx)*(tempt_x-seedx)+(tempt_y-seedy)*(tempt_y-seedy)+(tempt_z-seedz)*(tempt_z-seedz))
	UpdateHeap(tempt_x, tempt_y, tempt_z, MaxDist);
    }
  }


  /* ============================ */
  tempt_y=min(min_y+1,ydim1-1);
  tempt_x=min_x;
  tempt_z=min_z;
  if (heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)] == MaxVal) {
    dist = (tempt_x-min_seedx)*(tempt_x-min_seedx)+(tempt_y-min_seedy)*(tempt_y-min_seedy)+(tempt_z-min_seedz)*(tempt_z-min_seedz);
    if (dist <= threshold)
      InsertHeap(tempt_x, tempt_y, tempt_z, dist);
    else 
      boundary = 1;
  }
  else if (heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)] > -1) {
    neighbor = heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)];
    dt = min_heap->dist[neighbor];
    if (dt < MaxDist) {
      dist = (tempt_x-min_seedx)*(tempt_x-min_seedx) + (tempt_y-min_seedy)*(tempt_y-min_seedy) + (tempt_z-min_seedz)*(tempt_z-min_seedz);
      if (dist < dt)
	UpdateHeap(tempt_x, tempt_y, tempt_z, dist);
    }
    else {
      dist = (tempt_x-min_seedx)*(tempt_x-min_seedx) + (tempt_y-min_seedy)*(tempt_y-min_seedy) + (tempt_z-min_seedz)*(tempt_z-min_seedz);
      seed = min_heap->seed[neighbor];
      seedx = AllSeeds[seed].seedx;
      seedy = AllSeeds[seed].seedy;
      seedz = AllSeeds[seed].seedz;
      if (dist < (tempt_x-seedx)*(tempt_x-seedx)+(tempt_y-seedy)*(tempt_y-seedy)+(tempt_z-seedz)*(tempt_z-seedz))
	UpdateHeap(tempt_x, tempt_y, tempt_z, MaxDist);
    }
  }


  /* ============================ */
  tempt_z=max(min_z-1,0);
  tempt_x=min_x;
  tempt_y=min_y;
  if (heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)] == MaxVal) {
    dist = (tempt_x-min_seedx)*(tempt_x-min_seedx)+(tempt_y-min_seedy)*(tempt_y-min_seedy)+(tempt_z-min_seedz)*(tempt_z-min_seedz);
    if (dist <= threshold)
      InsertHeap(tempt_x, tempt_y, tempt_z, dist);
    else 
      boundary = 1;
  }
  else if (heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)] > -1) {
    neighbor = heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)];
    dt = min_heap->dist[neighbor];
    if (dt < MaxDist) {
      dist = (tempt_x-min_seedx)*(tempt_x-min_seedx) + (tempt_y-min_seedy)*(tempt_y-min_seedy) + (tempt_z-min_seedz)*(tempt_z-min_seedz);
      if (dist < dt)
	UpdateHeap(tempt_x, tempt_y, tempt_z, dist);
    }
    else {
      dist = (tempt_x-min_seedx)*(tempt_x-min_seedx) + (tempt_y-min_seedy)*(tempt_y-min_seedy) + (tempt_z-min_seedz)*(tempt_z-min_seedz);
      seed = min_heap->seed[neighbor];
      seedx = AllSeeds[seed].seedx;
      seedy = AllSeeds[seed].seedy;
      seedz = AllSeeds[seed].seedz;
      if (dist < (tempt_x-seedx)*(tempt_x-seedx)+(tempt_y-seedy)*(tempt_y-seedy)+(tempt_z-seedz)*(tempt_z-seedz))
	UpdateHeap(tempt_x, tempt_y, tempt_z, MaxDist);
    }
  }


  /* ============================ */
  tempt_z=min(min_z+1,zdim1-1);
  tempt_x=min_x;
  tempt_y=min_y;
  if (heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)] == MaxVal) {
    dist = (tempt_x-min_seedx)*(tempt_x-min_seedx)+(tempt_y-min_seedy)*(tempt_y-min_seedy)+(tempt_z-min_seedz)*(tempt_z-min_seedz);
    if (dist <= threshold)
      InsertHeap(tempt_x, tempt_y, tempt_z, dist);
    else 
      boundary = 1;
  }
  else if (heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)] > -1) {
    neighbor = heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)];
    dt = min_heap->dist[neighbor];
    if (dt < MaxDist) {
      dist = (tempt_x-min_seedx)*(tempt_x-min_seedx) + (tempt_y-min_seedy)*(tempt_y-min_seedy) + (tempt_z-min_seedz)*(tempt_z-min_seedz);
      if (dist < dt)
	UpdateHeap(tempt_x, tempt_y, tempt_z, dist);
    }
    else {
      dist = (tempt_x-min_seedx)*(tempt_x-min_seedx) + (tempt_y-min_seedy)*(tempt_y-min_seedy) + (tempt_z-min_seedz)*(tempt_z-min_seedz);
      seed = min_heap->seed[neighbor];
      seedx = AllSeeds[seed].seedx;
      seedy = AllSeeds[seed].seedy;
      seedz = AllSeeds[seed].seedz;
      if (dist < (tempt_x-seedx)*(tempt_x-seedx)+(tempt_y-seedy)*(tempt_y-seedy)+(tempt_z-seedz)*(tempt_z-seedz))
	UpdateHeap(tempt_x, tempt_y, tempt_z, MaxDist);
    }
  }


  if (boundary)
    InsertHeap(min_x,min_y,min_z, MaxDist);

}

 

FLTVECT FindSeed(float x, float y, float z, int index)
{
  double cx1,cy1,cz1;
  double cx2,cy2,cz2;
  double cx3,cy3,cz3;
  double dist,radius1,radius2,radius3;
  double cos_alpha;
  double ax,ay,az;
  double bx,by,bz;
  double dx,dy,dz;
  double hx,hy,hz;
  int atom1,atom2,atom3;
  int num,total;
  FLTVECT tmp;


  // AllSeeds[index].atom[0] is always >= 0
  atom1 = AllSeeds[index].atom[0];
  atom2 = AllSeeds[index].atom[1];
  atom3 = AllSeeds[index].atom[2];

  // contain only one atom
  if (atom2 < 0) {
    radius1 = atom_list[atom1].radius;
    cx1 = atom_list[atom1].x;
    cy1 = atom_list[atom1].y;
    cz1 = atom_list[atom1].z;
    dx = x-cx1;
    dy = y-cy1;
    dz = z-cz1;
    dist = sqrt(dx*dx + dy*dy + dz*dz);
    
    tmp.x = cx1 + radius1*dx/dist;
    tmp.y = cy1 + radius1*dy/dist;
    tmp.z = cz1 + radius1*dz/dist;
  }
  // contain only two atoms
  else if (atom3 < 0) {
    radius1 = atom_list[atom1].radius;
    cx1 = atom_list[atom1].x;
    cy1 = atom_list[atom1].y;
    cz1 = atom_list[atom1].z;
    radius2 = atom_list[atom2].radius;
    cx2 = atom_list[atom2].x;
    cy2 = atom_list[atom2].y;
    cz2 = atom_list[atom2].z;
    dist = sqrt((cx2-cx1)*(cx2-cx1)+(cy2-cy1)*(cy2-cy1)+(cz2-cz1)*(cz2-cz1));
    cos_alpha = (radius1*radius1+dist*dist-radius2*radius2)/(2.0*radius1*dist);
    
    ax = (cx2 - cx1)/dist;
    ay = (cy2 - cy1)/dist;
    az = (cz2 - cz1)/dist;
    bx = x-cx1;
    by = y-cy1;
    bz = z-cz1;
    dx = ay*bz - az*by;
    dy = az*bx - ax*bz;
    dz = ax*by - ay*bx;
    hx = dy*az - dz*ay;
    hy = dz*ax - dx*az;
    hz = dx*ay - dy*ax;
    dist = sqrt(hx*hx+hy*hy+hz*hz);
    hx /= dist;
    hy /= dist;
    hz /= dist;
    dist = radius1*sqrt(1.0-cos_alpha*cos_alpha);
    
    radius1 *= cos_alpha;
    tmp.x = radius1*ax+cx1 + dist*hx;
    tmp.y = radius1*ay+cy1 + dist*hy;
    tmp.z = radius1*az+cz1 + dist*hz;
  }
  // contain three or more atoms
  else {
    hx = 0;
    hy = 0;
    hz = 0;
    total = 0;
    num = 2;
    while (num < MaxAtom) {
      if (AllSeeds[index].atom[num] < 0)
	break;
      atom1 = AllSeeds[index].atom[num-2];
      atom2 = AllSeeds[index].atom[num-1];
      atom3 = AllSeeds[index].atom[num];
      
      radius1 = atom_list[atom1].radius;
      cx1 = atom_list[atom1].x;
      cy1 = atom_list[atom1].y;
      cz1 = atom_list[atom1].z;
      radius2 = atom_list[atom2].radius;
      cx2 = atom_list[atom2].x;
      cy2 = atom_list[atom2].y;
      cz2 = atom_list[atom2].z;
      radius3 = atom_list[atom3].radius;
      cx3 = atom_list[atom3].x;
      cy3 = atom_list[atom3].y;
      cz3 = atom_list[atom3].z;
      
      ax = x;
      ay = y;
      az = z;
      bx = ax;
      by = ay;
      bz = az;
      while (1) {
	// projection to the first atom sphere
	dx = ax-cx1;
	dy = ay-cy1;
	dz = az-cz1;
	dist = radius1/sqrt(dx*dx+dy*dy+dz*dz);
	ax = cx1 + dist*dx;
	ay = cy1 + dist*dy;
	az = cz1 + dist*dz;
	// projection to the second atom sphere
	dx = ax-cx2;
	dy = ay-cy2;
	dz = az-cz2;
	dist = radius2/sqrt(dx*dx+dy*dy+dz*dz);
	ax = cx2 + dist*dx;
	ay = cy2 + dist*dy;
	az = cz2 + dist*dz;
	// projection to the third atom sphere
	dx = ax-cx3;
	dy = ay-cy3;
	dz = az-cz3;
	dist = radius3/sqrt(dx*dx+dy*dy+dz*dz);
	ax = cx3 + dist*dx;
	ay = cy3 + dist*dy;
	az = cz3 + dist*dz;
	
	dist = sqrt((ax-bx)*(ax-bx)+(ay-by)*(ay-by)+(az-bz)*(az-bz));
	if (dist < 0.001) {
	  hx += ax;
	  hy += ay;
	  hz += az;
	  total++;
	  break;
	}
	else {
	  bx = ax;
	  by = ay;
	  bz = az;
	}
      }
      num++;
    }
    tmp.x = hx/(float)total;
    tmp.y = hy/(float)total;
    tmp.z = hz/(float)total;
  }
  
  return (tmp);
}

