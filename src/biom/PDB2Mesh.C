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
 * File:     PDB2Mesh.C    < ... >
 *
 * Author:   Zeyun Yu
 *
 * Purpose:  Convert PDB into 3D Volume using binary distance transform
 * ***************************************************************************
 */

#include <gamer/biom.h>
#include "gamercf.h"
#include <vector>
#include <cmath>

typedef struct {
  float x;       // vertex coordinate
  float y;
  float z;
  unsigned char neigh;    // first bit: +x; second bit: -x
                          // third bit: +y; fourth bit: -y
                          // fifth bit: +z; sixth  bit: -z
  unsigned short px;      // the corresponding index
  unsigned short py;
  unsigned short pz;
}MOL_VERTEX;

int *segment_index;
int *atom_index;
INT4VECT *quads;
MOL_VERTEX *vertex;
int vert_num;
int quad_num;
int xdim,ydim,zdim;
  

int CheckFaceCorner(float x, float y, float z);
void ExtractSES(MinHeapS*,SEEDS *, int *, int, int, int, int *, int, ATOM *, float);
void SetAtomIndex(int, int, int, int);
int ExtractSAS(int atom_num, ATOM *atom_list);
char CheckManifold(int i, int j, int k);
float GetAngle(int a, int b, int c);
void ReadPDB(char *filename, int *atom_num, ATOM **atomlist,
	     float min[3], float max[3]);

SurfaceMesh* SurfaceMesh_readPDB_molsurf(char* input_name)
{
  int i,j,k;
  int a,b,c,d;
  float orig[3],span[3];
  int dim[3];
  int m,n,l,num;
  double threshold;
  double nx,ny,nz;
  int xydim,xyzdim;
  clock_t begin, finish;
  MinHeapS* min_heap;
  SEEDS *AllSeeds;
  SurfaceMesh *surfmesh;
  int atom_num;
  ATOM *atom_list = NULL;
  float min[3], max[3];

  // Read in the PDB file
  ReadPDB(input_name, &atom_num, &atom_list, min, max);
  
  xdim = (int)(((max[0]-min[0])+1)*DIM_SCALE);
  ydim = (int)(((max[1]-min[1])+1)*DIM_SCALE);
  zdim = (int)(((max[2]-min[2])+1)*DIM_SCALE);
  xydim = xdim*ydim;
  xyzdim = xdim*ydim*zdim;
  //printf("dimension: %d X %d X %d\n",xdim,ydim,zdim);
  
  atom_index = (int*)malloc(sizeof(int)*xyzdim);
  segment_index = (int*)malloc(sizeof(int)*xyzdim);
  
  for (k=0; k<xyzdim; k++)
    atom_index[k] = 0;
  
  orig[0] = min[0];
  orig[1] = min[1];
  orig[2] = min[2];
  span[0] = (max[0] - min[0]) / (double)(xdim-1);
  span[1] = (max[1] - min[1]) / (double)(ydim-1);
  span[2] = (max[2] - min[2]) / (double)(zdim-1);
  dim[0] = xdim;
  dim[1] = ydim;
  dim[2] = zdim;

  for (m=0; m < atom_num; m++) {
    atom_list[m].x = (atom_list[m].x-orig[0])/span[0];
    atom_list[m].y = (atom_list[m].y-orig[1])/span[1];
    atom_list[m].z = (atom_list[m].z-orig[2])/span[2];
    atom_list[m].radius = (atom_list[m].radius+1.5)/((span[0]+span[1]+span[2])/3.0);
  }

  begin = clock();
  num = ExtractSAS(atom_num,atom_list);
  finish = clock();
  //printf("   Extract SAS voxels: CPU Time = %f seconds \n",(double)(finish-begin)/CLOCKS_PER_SEC);
  //printf("   Number of boundary voxels: %d\n\n",num);


  begin = clock();
  threshold = 1.5/((span[0]+span[1]+span[2])/3.0);
  min_heap=(MinHeapS*)malloc(sizeof(MinHeapS));
  min_heap->x = (unsigned short *)malloc(sizeof(unsigned short)*num*3);
  min_heap->y = (unsigned short *)malloc(sizeof(unsigned short)*num*3);
  min_heap->z = (unsigned short *)malloc(sizeof(unsigned short)*num*3);
  min_heap->seed = (int *)malloc(sizeof(int)*num*3);
  min_heap->dist = (float *)malloc(sizeof(float)*num*3);
  AllSeeds = (SEEDS *)malloc(sizeof(SEEDS)*num);
  ExtractSES(min_heap, AllSeeds, segment_index, xdim,ydim,zdim,atom_index,atom_num,atom_list, threshold*threshold);
  finish = clock();
  //printf("   Extract SES voxels: CPU Time = %f seconds \n\n",(double)(finish-begin)/CLOCKS_PER_SEC);
  

  
  // detect and fix non-manifolds !
  while (1) {
    b = 0;
    for (k=0; k<zdim; k++)
      for (j=0; j<ydim; j++)
	for (i=0; i<xdim; i++) {
	  if (segment_index[IndexVect(i,j,k)] == MaxVal) { 
	    if (!CheckManifold(i, j, k)) { // non-manifold occurs 
	      segment_index[IndexVect(i,j,k)] = 0;
	      min_heap->x[min_heap->size] = i;
	      min_heap->y[min_heap->size] = j;
	      min_heap->z[min_heap->size] = k;
	      min_heap->size++;
	      b++;
	    }
	  }
	}
    //printf("non-manifold = %d \n",b);
    if (b == 0)
      break;
  }
  

	    
  // generate the surface mesh
  vertex = (MOL_VERTEX *)malloc(sizeof(MOL_VERTEX)*num*8);
  quads = (INT4VECT *)malloc(sizeof(INT4VECT)*num*6);
  for (k=0; k<num*8; k++) {
    vertex[k].neigh = 0;
  }
  for (k=0; k<xyzdim; k++)
    atom_index[k] = -1;
  begin = clock();
  vert_num = 0;
  quad_num = 0;
  for(num = 0; num < min_heap->size; num++) {
    i = min_heap->x[num];
    j = min_heap->y[num];
    k = min_heap->z[num];
    
    // back face
    if (segment_index[IndexVect(i-1,j,k)] == MaxVal) {
      a = CheckFaceCorner(i-0.5, j-0.5, k-0.5);
      vertex[a].neigh |= 40; // +y and +z
      b = CheckFaceCorner(i-0.5, j-0.5, k+0.5);
      vertex[b].neigh |= 24; // +y and -z
      c = CheckFaceCorner(i-0.5, j+0.5, k+0.5);
      vertex[c].neigh |= 20; // -y and -z
      d = CheckFaceCorner(i-0.5, j+0.5, k-0.5);
      vertex[d].neigh |= 36; // -y and +z
	        
      quads[quad_num].a = a;
      quads[quad_num].b = b;
      quads[quad_num].c = c;
      quads[quad_num].d = d;
      quad_num++;
    }
    
    // front face
    if (segment_index[IndexVect(i+1,j,k)] == MaxVal) {
      a = CheckFaceCorner(i+0.5, j-0.5, k-0.5);
      vertex[a].neigh |= 40; // +y and +z
      b = CheckFaceCorner(i+0.5, j+0.5, k-0.5);
      vertex[b].neigh |= 36; // -y and +z
      c = CheckFaceCorner(i+0.5, j+0.5, k+0.5);
      vertex[c].neigh |= 20; // -y and -z
      d = CheckFaceCorner(i+0.5, j-0.5, k+0.5);
      vertex[d].neigh |= 24; // +y and -z
      
      quads[quad_num].a = a;
      quads[quad_num].b = b;
      quads[quad_num].c = c;
      quads[quad_num].d = d;
      quad_num++;
    }
    
    // left face
    if (segment_index[IndexVect(i,j-1,k)] == MaxVal) {
      a = CheckFaceCorner(i+0.5, j-0.5, k-0.5);
      vertex[a].neigh |= 33; // -x and +z
      b = CheckFaceCorner(i+0.5, j-0.5, k+0.5);
      vertex[b].neigh |= 17; // -x and -z
      c = CheckFaceCorner(i-0.5, j-0.5, k+0.5);
      vertex[c].neigh |= 18; // +x and -z
      d = CheckFaceCorner(i-0.5, j-0.5, k-0.5);
      vertex[d].neigh |= 34; // +x and +z

      quads[quad_num].a = a;
      quads[quad_num].b = b;
      quads[quad_num].c = c;
      quads[quad_num].d = d;
      quad_num++;
    }
    
    // right face
    if (segment_index[IndexVect(i,j+1,k)] == MaxVal) {
      a = CheckFaceCorner(i+0.5, j+0.5, k-0.5);
      vertex[a].neigh |= 33; // -x and +z
      b = CheckFaceCorner(i-0.5, j+0.5, k-0.5);
      vertex[b].neigh |= 34; // +x and +z
      c = CheckFaceCorner(i-0.5, j+0.5, k+0.5);
      vertex[c].neigh |= 18; // +x and -z
      d = CheckFaceCorner(i+0.5, j+0.5, k+0.5);
      vertex[d].neigh |= 17; // -x and -z

      quads[quad_num].a = a;
      quads[quad_num].b = b;
      quads[quad_num].c = c;
      quads[quad_num].d = d;
      quad_num++;
    }
    
    // bottom face
    if (segment_index[IndexVect(i,j,k-1)] == MaxVal) {
      a = CheckFaceCorner(i+0.5, j-0.5, k-0.5);
      vertex[a].neigh |= 9; // -x and +y
      b = CheckFaceCorner(i-0.5, j-0.5, k-0.5);
      vertex[b].neigh |= 10; // +x and +y
      c = CheckFaceCorner(i-0.5, j+0.5, k-0.5);
      vertex[c].neigh |= 6; // +x and -y
      d = CheckFaceCorner(i+0.5, j+0.5, k-0.5);
      vertex[d].neigh |= 5; // -x and -y
      
      quads[quad_num].a = a;
      quads[quad_num].b = b;
      quads[quad_num].c = c;
      quads[quad_num].d = d;
      quad_num++;
    }
    
    // top face
    if (segment_index[IndexVect(i,j,k+1)] == MaxVal) {
      a = CheckFaceCorner(i+0.5, j-0.5, k+0.5);
      vertex[a].neigh |= 9; // -x and +y
      b = CheckFaceCorner(i+0.5, j+0.5, k+0.5);
      vertex[b].neigh |= 5; // -x and -y
      c = CheckFaceCorner(i-0.5, j+0.5, k+0.5);
      vertex[c].neigh |= 6; // +x and -y
      d = CheckFaceCorner(i-0.5, j-0.5, k+0.5);
      vertex[d].neigh |= 10; // +x and +y

      quads[quad_num].a = a;
      quads[quad_num].b = b;
      quads[quad_num].c = c;
      quads[quad_num].d = d;
      quad_num++;
    }
  }
  finish = clock();
  //printf("   Generate quad meshes: CPU Time = %f seconds \n",(double)(finish-begin)/CLOCKS_PER_SEC);
  //printf("   vert-num : %d -- quad-num: %d \n\n",vert_num,quad_num);
	    
  
  
  // Smooth the mesh 
  begin = clock();
  unsigned char neighbor;
  for (num = 0; num < 3; num++) {
  for (n = 0; n < vert_num; n++) {
    nx = 0;
    ny = 0;
    nz = 0;
    m = 0;
    neighbor = vertex[n].neigh;

    i = vertex[n].px;
    j = vertex[n].py;
    k = vertex[n].pz;
    if (neighbor & 1) {
      m++;
      l = atom_index[IndexVect(i-1,j,k)];
      nx += vertex[l].x;
      ny += vertex[l].y;
      nz += vertex[l].z;
    }
    if (neighbor & 2) {
      m++;
      l = atom_index[IndexVect(i+1,j,k)];
      nx += vertex[l].x;
      ny += vertex[l].y;
      nz += vertex[l].z;
    }
    if (neighbor & 4) {
      m++;
      l = atom_index[IndexVect(i,j-1,k)];
      nx += vertex[l].x;
      ny += vertex[l].y;
      nz += vertex[l].z;
    }
    if (neighbor & 8) {
      m++;
      l = atom_index[IndexVect(i,j+1,k)];
      nx += vertex[l].x;
      ny += vertex[l].y;
      nz += vertex[l].z;
    }
    if (neighbor & 16) {
      m++;
      l = atom_index[IndexVect(i,j,k-1)];
      nx += vertex[l].x;
      ny += vertex[l].y;
      nz += vertex[l].z;
    }
    if (neighbor & 32) {
      m++;
      l = atom_index[IndexVect(i,j,k+1)];
      nx += vertex[l].x;
      ny += vertex[l].y;
      nz += vertex[l].z;
    }
    
    // update the position
    vertex[n].x = nx/(float)m;
    vertex[n].y = ny/(float)m;
    vertex[n].z = nz/(float)m;
  }
  }
  finish = clock();
  //printf("   Smooth the quad meshes: CPU Time = %f seconds \n\n",(double)(finish-begin)/CLOCKS_PER_SEC);
  
  // Allocate memory
  surfmesh = SurfaceMesh_ctor(vert_num, quad_num*2);

  // write vertices
  for (i = 0; i < surfmesh->num_vertices; i++) {
    surfmesh->vertex[i].x = vertex[i].x*span[0]+orig[0];
    surfmesh->vertex[i].y = vertex[i].y*span[1]+orig[1];
    surfmesh->vertex[i].z = vertex[i].z*span[2]+orig[2];
  }
  // write triangles
  float angle, angle1, angle2;
  for (i = 0; i < quad_num; i++) {
    a = quads[i].a;
    b = quads[i].b;
    c = quads[i].c;
    d = quads[i].d;
    
    angle1 = -999.0;
    angle2 = -999.0;
    angle = GetAngle(a,b,c);
    if (angle > angle1)
      angle1 = angle;
    angle = GetAngle(a,d,c);
    if (angle > angle1)
      angle1 = angle;
    angle = GetAngle(c,a,b);
    if (angle > angle1)
      angle1 = angle;
    angle = GetAngle(c,a,d);
    if (angle > angle1)
      angle1 = angle;
    
    angle = GetAngle(b,a,d);
    if (angle > angle2)
      angle2 = angle;
    angle = GetAngle(b,c,d);
    if (angle > angle2)
      angle2 = angle;
    angle = GetAngle(d,a,b);
    if (angle > angle2)
      angle2 = angle;
    angle = GetAngle(d,b,c);
    if (angle > angle2)
      angle2 = angle;
    
    if (angle1 <= angle2) {
      surfmesh->face[2*i].a = a;
      surfmesh->face[2*i].b = b;
      surfmesh->face[2*i].c = c;
      surfmesh->face[2*i+1].a = a;
      surfmesh->face[2*i+1].b = c;
      surfmesh->face[2*i+1].c = d;
    }
    else {
      surfmesh->face[2*i].a = a;
      surfmesh->face[2*i].b = b;
      surfmesh->face[2*i].c = d;
      surfmesh->face[2*i+1].a = b;
      surfmesh->face[2*i+1].b = c;
      surfmesh->face[2*i+1].c = d;
    }
  }

  /* write to disk 
  FILE *fout;
  if ((fout=fopen("output1.off", "wb"))==NULL){
    printf("write error...\n");
    exit(0);
  };
  
  fprintf(fout, "OFF\n");
  fprintf(fout, "%d %d %d\n",vert_num,quad_num*2,vert_num+quad_num*2-2);
  for (i = 0; i < vert_num; i++)
    fprintf(fout, "%f %f %f \n",surfmesh->vertex[i].x,surfmesh->vertex[i].y,surfmesh->vertex[i].z);

  for (i = 0; i < quad_num*2; i++) {
    fprintf(fout, "3 %d %d %d\n",surfmesh->face[i].a,surfmesh->face[i].b,surfmesh->face[i].c);
  }
    fclose(fout);
  */

  free(AllSeeds);
  free(segment_index);
  free(atom_index);
  free(atom_list);
  free(min_heap->x);
  free(min_heap->y);
  free(min_heap->z);
  free(min_heap->seed);
  free(min_heap->dist);
  free(min_heap);
  free(vertex);
  free(quads);

  // Split multiple connectedsurfaces
  SurfaceMesh_splitMultipleConnectedSurfaces(surfmesh);

  // Flip normals so they now points outwards
  SurfaceMesh_flipNormals(surfmesh);
  
  return surfmesh;
}

float GetAngle(int a, int b, int c)
{
  float ax,ay,az;
  float bx,by,bz;
  float dist;

  ax = vertex[b].x - vertex[a].x;
  ay = vertex[b].y - vertex[a].y;
  az = vertex[b].z - vertex[a].z;
  dist = sqrt(ax*ax+ay*ay+az*az);
  if (dist > 0) {
    ax /= dist;
    ay /= dist;
    az /= dist;
  }
  bx = vertex[c].x - vertex[a].x;
  by = vertex[c].y - vertex[a].y;
  bz = vertex[c].z - vertex[a].z;
  dist = sqrt(bx*bx+by*by+bz*bz);
  if (dist > 0) {
    bx /= dist;
    by /= dist;
    bz /= dist;
  }

  return (ax*bx+ay*by+az*bz);
}




char CheckManifold(int i, int j, int k)
{
  char manifold, nonmanifold;
  int m,n,l;

  manifold = 1;
  for (l=k-1; l<=k+1; l++) {
    if (!manifold)
      break;
    for (n=j-1; n<=j+1; n++) {
      if (!manifold)
	break;
      for (m=i-1; m<=i+1; m++) {
	if (m!=i || n!=j || l!=k) {
	  if (segment_index[IndexVect(m,n,l)] == MaxVal) { 
	    nonmanifold = 1;
	    if (m!=i && segment_index[IndexVect(m,j,k)] == MaxVal)
	      nonmanifold = 0;
	    if (n!=j && segment_index[IndexVect(i,n,k)] == MaxVal)
	      nonmanifold = 0;
	    if (l!=k && segment_index[IndexVect(i,j,l)] == MaxVal)
	      nonmanifold = 0;
	    
	    if (nonmanifold) 
	      manifold = 0;
	  }
	}
	if (!manifold)
	  break;
      }
    }
  }
  
  return(manifold);
}

  

int CheckFaceCorner(float x, float y, float z)
{
  int m,n,l;
  int a;

  m = (int)x;
  n = (int)y;
  l = (int)z;
  if (atom_index[IndexVect(m,n,l)] < 0) {
    vertex[vert_num].x = x;
    vertex[vert_num].y = y;
    vertex[vert_num].z = z;
    vertex[vert_num].px = m;
    vertex[vert_num].py = n;
    vertex[vert_num].pz = l;
    atom_index[IndexVect(m,n,l)] = vert_num;
    a = vert_num;
    vert_num++;
  }
  else {
    a = atom_index[IndexVect(m,n,l)];
  }

  return (a);
}



FLT2VECT FindIntersection(int n, int m, int j, int k, ATOM *atom_list)
{
  FLT2VECT intersect;
  int i;
  char reorder;
  double cx1,cy1,cz1;
  double cx2,cy2,cz2;
  double dist,radius1,radius2;
  double cos_alpha;
  double ax,ay,az;
  double cx,cy,cz;


  n--;
  m--;
  reorder = 0;
  // always make the first atom on the left
  // and the second on the right
  if (atom_list[n].x > atom_list[m].x) {
    i = n;
    n = m;
    m = i;
    reorder = 1;
  }

  radius1 = atom_list[n].radius;
  cx1 = atom_list[n].x;
  cy1 = atom_list[n].y;
  cz1 = atom_list[n].z;
  radius2 = atom_list[m].radius;
  cx2 = atom_list[m].x;
  cy2 = atom_list[m].y;
  cz2 = atom_list[m].z;
  
  if (cx1 == cx2) {
    if ((radius1*radius1-(k-cz1)*(k-cz1)-(j-cy1)*(j-cy1)) >=
        (radius2*radius2-(k-cz2)*(k-cz2)-(j-cy2)*(j-cy2))) {
      intersect.x = (float)MaxVal;
      intersect.y = -(float)MaxVal;
    }
    else {
      intersect.x = -(float)MaxVal;
      intersect.y = (float)MaxVal;
    }
  }
  else {
    dist = sqrt((cx2-cx1)*(cx2-cx1)+(cy2-cy1)*(cy2-cy1)+(cz2-cz1)*(cz2-cz1));
    cos_alpha = (radius1*radius1+dist*dist-radius2*radius2)/(2.0*radius1*dist);
    ax = (cx2 - cx1)/dist;
    ay = (cy2 - cy1)/dist;
    az = (cz2 - cz1)/dist;
    dist = radius1 * cos_alpha;
    cx = dist*ax+cx1;
    cy = dist*ay+cy1;
    cz = dist*az+cz1;
    dist = cx*ax + (cy-j)*ay + (cz-k)*az;

    if (reorder == 0) {
      intersect.x = (float)(dist/ax);
      intersect.y = (float)MaxVal;
    }
    else {
      intersect.x = -(float)MaxVal;
      intersect.y = (float)(dist/ax);
    }
  }

  return (intersect);
}



int ExtractSAS(int atom_num, ATOM *atom_list)
{
  int i,j,k;
  int m,n,l;
  int dim[3], c[3];
  int amax[3],amin[3];
  float radius;
  float x,y,z;
  FLT2VECT intersect;


  dim[0] = xdim;
  dim[1] = ydim;
  dim[2] = zdim;


  for (m=0; m < atom_num; m++) {
    
    radius = atom_list[m].radius;
    
    // compute the dataset coordinates of the atom's center
    c[0] = (int)(atom_list[m].x + 0.5);
    c[1] = (int)(atom_list[m].y + 0.5);
    c[2] = (int)(atom_list[m].z + 0.5);
    
    // then compute the bounding box of the atom
    for (j=0; j < 3; j++) {
      int tmp;
      tmp = (int)(c[j] - radius - 1);
      tmp = (tmp < 0) ? 0 : tmp;
      amin[j] = tmp;
      tmp = (int)(c[j] + radius + 1);
      tmp = (tmp > (dim[j]-1)) ? (dim[j]-1) : tmp;
      amax[j] = tmp;
    }
    
    // begin blurring kernel
    radius = radius*radius;
    for(k = amin[2]; k <= amax[2]; k++) {
      for(j = amin[1]; j <= amax[1]; j++) {
        for(i = amin[0]; i <= amax[0]; i++) {
          x = i;
          y = j;
          z = k;
          if ((x-atom_list[m].x)*(x-atom_list[m].x)+
              (y-atom_list[m].y)*(y-atom_list[m].y)+
              (z-atom_list[m].z)*(z-atom_list[m].z) <= radius) {
	    if (atom_index[IndexVect(i,j,k)] > 0) {
	      intersect = FindIntersection(atom_index[IndexVect(i,j,k)], m+1, j, k, atom_list);
	      if (i >= intersect.x && i <= intersect.y)
		atom_index[IndexVect(i,j,k)] = m+1;
	    }
	    else
	      atom_index[IndexVect(i,j,k)] = m+1;
	  }
	}
      }
    }
    
  }
  
  
  // find voxels on the border
  int total = 0;
  for (l=1; l<zdim-1; l++) 
    for (n=1; n<ydim-1; n++)
      for (m=1; m<xdim-1; m++) {
        if (atom_index[IndexVect(m,n,l)]) {
          int count = 0;
          for (k=std::max(l-1,0); k<=std::min(l+1,zdim-1); k++) 
            for (j=std::max(n-1,0); j<=std::min(n+1,ydim-1); j++)
              for (i=std::max(m-1,0); i<=std::min(m+1,xdim-1); i++) {
                if (((i==m&&j==n) || (i==m&&k==l) || (k==l&&j==n) ) &&
		    atom_index[IndexVect(i,j,k)] == 0)
                  count = 1;
              }
          if (count) {
	    atom_index[IndexVect(m,n,l)] = -atom_index[IndexVect(m,n,l)];
            total++;
	  }
        }
      }
  
  
  return(total);
}
