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
 * File:     PDB2Volume.C    < ... >
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Convert a PDB into a 3D Volume
 *           
 * Source:   Part of this file was adapted from the PDBParser developed in 
 *           Chandrajit Bajaj's group at The University of Texas at Austin.
 * ***************************************************************************
 */


#include "PDB2Volume.h"

#define EPSILON		 1.0e-3f
#define MAX_STRING       256 


void getMinMax(ATOM *, int, float*, float*);
float evalDensity(ATOM *,int, float*, double);
void blurAtoms(ATOM *,int,float*,float*,float *,int*);
void write_rawiv_float(FILE*,float*,int*,float*,float*);
float PDB2Volume(char*, float**, int*, int*, int*, float[3], float[3], ATOM**, 
		 int*, char);

/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_readPDB_gauss    < ... >
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Convert a PDB/PQR/XYZR formats into a SurfaceMesh using Gaussian 
 *           scalar function
 *
 * Notes:    Each atom is treated as a Gaussian function defined by 
 *           - center and radius: (available in the PDB/PQR files)
 *           - blobbyness: the decay rate of the Gaussian function
 * ***************************************************************************
 */
SurfaceMesh* SurfaceMesh_readPDB_gauss(char* filename, float blobbyness, 
				       float iso_value)
{
  float max_density;
  float *dataset;
  float data_iso_val;
  time_t t1,t2;
  int xdim,ydim,zdim;
  int atom_num = 0;
  ATOM *atom_list = NULL;
  char IsXYZR;
  float min[3],max[3],span[3];
  SurfaceMesh *surfmesh;
  SPNT *holelist;
  int i,j,k;

  printf("\nbegin blurring PDB/PQR coordinates ... \n");
  (void)time(&t1);
  max_density = PDB2Volume(filename, &dataset, &xdim, &ydim, &zdim, min,max,
			   &atom_list, &atom_num, IsXYZR);
  (void)time(&t2);
  printf("time to generate volume from PDB: %d seconds. \n\n",(int)(t2-t1));
  span[0] = (max[0]-min[0])/(float)(xdim-1);
  span[1] = (max[1]-min[1])/(float)(ydim-1);
  span[2] = (max[2]-min[2])/(float)(zdim-1);
  
  printf("begin extracting isosurfaces ... \n");
  (void)time(&t1);
  data_iso_val = 0.44*max_density;
  if (data_iso_val < iso_value)
    iso_value = data_iso_val;
  printf("isovalue: %f \n", iso_value);
  surfmesh = SurfaceMesh_marchingCube(xdim, ydim, zdim, dataset, 
				      iso_value, &holelist);
  (void)time(&t2);
  printf("vertices: %d, faces: %d\n",surfmesh->num_vertices, surfmesh->num_faces);
  printf("time to extract isosurface: %d seconds. \n\n",(int)(t2-t1));
  free(dataset);
  
  // convert from pixel to angstrom
  for (j=0; j<surfmesh->num_vertices; j++) {
    surfmesh->vertex[j].x = surfmesh->vertex[j].x*span[0]+min[0];
    surfmesh->vertex[j].y = surfmesh->vertex[j].y*span[1]+min[1];
    surfmesh->vertex[j].z = surfmesh->vertex[j].z*span[2]+min[2];
  }

  // Flip normals so they now points outwards
  SurfaceMesh_flipNormals(surfmesh);
  
  // Return generated SurfaceMesh
  return surfmesh;
}

/*
 * ***************************************************************************
 * Routine:  PDB2Volume    < ... >
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Convert a PDB/PQR/XYZR formats into a 3D Volume
 *
 * Notes:    Each atom is treated as a Gaussian function defined by 
 *           - center and radius: (available in the PDB/PQR files)
 *           - blobbyness: the decay rate of the Gaussian function
 * ***************************************************************************
 */
float PDB2Volume(char *filename, float **data, int *xd, int *yd, int *zd,
		 float p_min[3], float p_max[3],ATOM **atomlist, int *atom_num, char IsXYZR)
{
  int m,n,k;
  char line[MAX_STRING];
  ATOM *atom_list;
  char string[8];
  float min[3],max[3];
  int dim[3];
  float min_dimension;
  float *dataset;
  FILE *fp,*fout;
  PDBelementInformation eInfo;
  char IsPQR,IsPDB;
  float x,y,z,radius;
  char file_name[256];


  if (IsXYZR) {
    if ((fp=fopen(filename, "r"))==NULL){
      printf("read error...\n");
      exit(0);
    };
    fscanf(fp,"%d\n",&m);
    atom_list = (ATOM*)malloc(sizeof(ATOM)*m);
    for (n = 0; n < m; n++) {
      fscanf(fp,"%f %f %f %f\n",&x,&y,&z,&radius);
      
      atom_list[n].x = x;
      atom_list[n].y = y;
      atom_list[n].z = z;
      atom_list[n].radius = radius;
      
    }
    
    *atom_num = m;
    
    printf("number of atoms: %d\n",m);
  }
  else {
    IsPQR = 0;
    IsPDB = 0;
    for(m = 0; m<256; m++) {
      if (filename[m+3] == '\0')
	break;
      else if (filename[m] == '.' && 
	       (filename[m+1] == 'P' || filename[m+1] == 'p') && 
	       (filename[m+2] == 'Q' || filename[m+2] == 'q') &&
	       (filename[m+3] == 'R' || filename[m+3] == 'r')/* &&
								file_name[m+4] == '\0'*/) {
	IsPQR = 1;
	break;
      }
      else if (filename[m] == '.' &&
	       (filename[m+1] == 'P' || filename[m+1] == 'p') &&
	       (filename[m+2] == 'D' || filename[m+2] == 'd') &&
	       (filename[m+3] == 'B' || filename[m+3] == 'b') /*&&
								file_name[m+4] == '\0'*/) {
	IsPDB = 1;
	break;
      }
    }
    if (IsPQR == 0 && IsPDB == 0) {
      printf("Input file name must be ending with PDB/PQR/XYZR/RAWIV/OFF...\n");
      exit(0);
    }

    if (IsPQR) {
      sprintf(file_name, "%s.xyzr", filename);
      if ((fout=fopen(file_name, "wb"))==NULL){
	printf("write error...\n");
	exit(0); 
      } 
    }
    
    if ((fp=fopen(filename, "r"))==NULL){
      printf("read error...\n");
      exit(0); 
    };
    m = 0;
    while (fgets(line,MAX_STRING,fp) != NULL) {
      if (line[0]=='A' && line[1]=='T' && line[2]=='O' && line[3]=='M') 
	m++;
    }
    printf("number of atoms: %d \n",m);
    fclose(fp);
    
    *atom_num = m;
    atom_list = (ATOM*)malloc(sizeof(ATOM)*m);
    if ((fp=fopen(filename, "r"))==NULL){
      printf("read error...\n");
      exit(0); 
    };
    m = 0;
    while (fgets(line,MAX_STRING,fp) != NULL) {
      if (line[0]=='A' && line[1]=='T' && line[2]=='O' && line[3]=='M') {
	/* more general format, could be used for pqr format */
	k = 30;
	while (line[k] == ' ') k++;
	n = 0;
	while (line[k] != ' ') {
	  string[n] = line[k];
	  n++;
	  k++;
	  if (line[k] == '-') 
	    break;
	}
	string[n] = '\0';
	atom_list[m].x = atof(string);
	while (line[k] == ' ') k++;
	n = 0;
	while (line[k] != ' ') {
	  string[n] = line[k];
	  n++;
	  k++;
	  if (line[k] == '-') 
	    break;
	}
	string[n] = '\0';
	atom_list[m].y = atof(string);
	while (line[k] == ' ') k++;
	n = 0;
	while (line[k] != ' ') {
	  string[n] = line[k];
	  n++;
	  k++;
	  if (line[k] == '-') 
	    break;
	}
	string[n] = '\0';
	atom_list[m].z = atof(string);
	
	if (IsPDB) {
	  atom_list[m].radius = 1.0f; // default radius
	  for( n=0; n<MAX_BIOCHEM_ELEMENTS; n++ ) {
	    eInfo = PDBelementTable[n];
	    if( (eInfo.atomName[0] == line[12]) &&  
		(eInfo.atomName[1] == line[13]) &&
		(eInfo.atomName[2] == line[14]) &&
		(eInfo.atomName[3] == line[15]) &&
		(eInfo.residueName[0] == line[17]) &&
		(eInfo.residueName[1] == line[18]) &&
		(eInfo.residueName[2] == line[19]) ) {
	      atom_list[m].radius = eInfo.radius;
	      break;
	    }
	  }
	}
	else if (IsPQR) {
	  while (line[k] == ' ') k++;
	  n = 0;
	  while (line[k] != ' ') {
	    string[n] = line[k];
	    n++;
	    k++;
	    if (line[k] == '-') 
	      break;
	  }
	  while (line[k] == ' ') k++;
	  n = 0;
	  while (line[k] != ' ' && line[k] != '\n' && line[k] != '\0') {
	    string[n] = line[k];
	    n++;
	    k++;
	  } 
	  string[n] = '\0';
	  atom_list[m].radius = atof(string);
	  //atom_list[m].radius += 1.0;
	 
	  if (atom_list[m].radius < 1.0)
	    atom_list[m].radius = 1.0;
	  
	  fprintf(fout, "%f %f %f %f\n",atom_list[m].x,atom_list[m].y,atom_list[m].z,atom_list[m].radius);
	}
	else {
	  printf("Input file name must be ending with PDB or PQR..\n");
	  exit(0);
	}

	m++;
      }
    }	
    fclose(fp);
  }
  if (IsPQR)
    fclose(fout);


  min[0] = min[1] = min[2] = 0.;
  max[0] = max[1] = max[2] = 0.;
  getMinMax(atom_list, m, min, max);
  p_min[0] = min[0];
  p_min[1] = min[1];
  p_min[2] = min[2];
  p_max[0] = max[0];
  p_max[1] = max[1];
  p_max[2] = max[2];

  min_dimension = min((max[0]-min[0]),min((max[1]-min[1]),(max[2]-min[2])));
  if (min_dimension < 64.0f) {
    min_dimension = 64.0f/min_dimension;
    dim[0] = (int)((max[0]-min[0])*min_dimension)+1;
    dim[1] = (int)((max[1]-min[1])*min_dimension)+1;
    dim[2] = (int)((max[2]-min[2])*min_dimension)+1;
  }
  else {
    dim[0] = (int)(max[0]-min[0])+1;
    dim[1] = (int)(max[1]-min[1])+1;
    dim[2] = (int)(max[2]-min[2])+1;
  }
  dim[0] = (int)(dim[0]*DIM_SCALE);
  dim[1] = (int)(dim[1]*DIM_SCALE);
  dim[2] = (int)(dim[2]*DIM_SCALE);

  printf("dimension: %d X %d X %d\n",dim[0],dim[1],dim[2]);

  dataset = (float*)malloc(sizeof(float)*dim[0]*dim[1]*dim[2]);
  blurAtoms(atom_list, m, min, max, dataset, dim);
  printf("min[3]: %f %f %f \n",min[0],min[1],min[2]);
  printf("max[3]: %f %f %f \n",max[0],max[1],max[2]);
  printf("span[3]: %f %f %f \n",(max[0] - min[0]) / (float)(dim[0]-1),
	 (max[1]-min[1])/(float)(dim[1]-1),(max[2]-min[2])/(float)(dim[2]-1));

  float minval,maxval;
  minval = 999999.0;
  maxval = -999999.0;
  for (k = 0; k < dim[2]; k++)
    for (n = 0; n < dim[1]; n++)
      for (m = 0; m < dim[0]; m++) {
	if (dataset[((k)*dim[1]*dim[0]+(n)*dim[0]+(m))] < minval)
	  minval = dataset[((k)*dim[1]*dim[0]+(n)*dim[0]+(m))];
	if (dataset[((k)*dim[1]*dim[0]+(n)*dim[0]+(m))] > maxval)
	  maxval = dataset[((k)*dim[1]*dim[0]+(n)*dim[0]+(m))];
      }
  printf("min_density: %f   max_density: %f \n",minval,maxval);


  /* write the volume to disk 
  if ((fp=fopen("test.rawiv", "w"))==NULL){
    printf("write error...\n");
    exit(0); 
  };
  write_rawiv_float(fp,dataset,dim,min,max);
  exit(0);
  */


  *data = dataset;
  *xd = dim[0];
  *yd = dim[1];
  *zd = dim[2];
  *atomlist = atom_list;
  return(maxval);
}



/*
 * ***************************************************************************
 * Routine:  blurAtoms    < ... >
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Blur all atoms with Gaussian functions
 *           Each atom spreads in a region of certain radius 
 * ***************************************************************************
 */
void blurAtoms(ATOM *atom_list, int size, float min[3], float max[3],
	       float *dataset, int dim[3])
{
  float orig[3],span[3];
  int i,j,k;
  int m;
  int xdim,ydim,zdim;
  int amax[3],amin[3];
  double c[3], maxRad;
  
  
  xdim = dim[0];
  ydim = dim[1];
  zdim = dim[2];
  for (k=0; k<zdim; k++) 
    for (j=0; j<ydim; j++)
      for (i=0; i<xdim; i++) 
	dataset[IndexVect(i,j,k)] = 0;
  
  orig[0] = min[0];
  orig[1] = min[1];
  orig[2] = min[2];
  span[0] = (max[0] - min[0]) / (float)(xdim-1);
  span[1] = (max[1] - min[1]) / (float)(ydim-1);
  span[2] = (max[2] - min[2]) / (float)(zdim-1);
  
  for (m=0; m < size; m++) {
    
    maxRad = atom_list[m].radius * sqrt(1.0 + log(EPSILON) / (2.0*BLOBBYNESS) );

    // compute the dataset coordinates of the atom's center
    c[0] = (atom_list[m].x-orig[0])/span[0];
    c[0] = ((c[0]-floor(c[0])) >= 0.5) ? ceil(c[0]) : floor(c[0]);
    c[1] = (atom_list[m].y-orig[1])/span[1];
    c[1] = ((c[1]-floor(c[1])) >= 0.5) ? ceil(c[1]) : floor(c[1]);
    c[2] = (atom_list[m].z-orig[2])/span[2];
    c[2] = ((c[2]-floor(c[2])) >= 0.5) ? ceil(c[2]) : floor(c[2]);
    
    // then compute the bounding box of the atom (maxRad^3)
    for (j=0; j < 3; j++) {
      int tmp;
      tmp = (int)(c[j] - (maxRad / span[j]) - 1);
      tmp = (tmp < 0) ? 0 : tmp;
      amin[j] = tmp;
      tmp = (int)(c[j] + (maxRad / span[j]) + 1);
      tmp = (tmp > (dim[j]-1)) ? (dim[j]-1) : tmp;
      amax[j] = tmp;
    }

    // begin blurring kernel
    for(k = amin[2]; k <= amax[2]; k++) {
      for(j = amin[1]; j <= amax[1]; j++) {
	for(i = amin[0]; i <= amax[0]; i++) {
	  float pnt[3], density;
	  
	  pnt[0] = orig[0] + i*span[0];
	  pnt[1] = orig[1] + j*span[1];
	  pnt[2] = orig[2] + k*span[2];

	  density = evalDensity(atom_list,m, pnt, maxRad);
	  dataset[IndexVect(i,j,k)] += density;
	}
      }
    }

    if (((m+1) % 20) == 0 || (m+1) == size) {
      printf("%2.2f%% done (%08d)\r", 100.0*(m+1)/(float)size, m+1);
      fflush(stdout);
      }
  }
  printf("\n"); fflush(stdout);
}




/*
 * ***************************************************************************
 * Routine:  evalDensity    < ... >
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Calculate the Gaussian function value using the blobbyness
 * ***************************************************************************
 */
float evalDensity(ATOM *atom_list,int index, float pnt[3], double maxRadius)
{
  double expval;

  double r = (atom_list[index].x-pnt[0])*(atom_list[index].x-pnt[0]) +
    (atom_list[index].y-pnt[1])*(atom_list[index].y-pnt[1]) +
    (atom_list[index].z-pnt[2])*(atom_list[index].z-pnt[2]);
  double r0 = atom_list[index].radius; r0 *= r0;
  //expval = BLOBBYNESS*(r/r0 - 1.0);
  expval = BLOBBYNESS*(r-r0);
  
  // truncated gaussian
  if (sqrt(r) > maxRadius)
    return 0.0;
  
  return (float)(exp(expval));
}




/*
 * ***************************************************************************
 * Routine:  getMinMax    < ... >
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Calculate the minimum and maximum range of the blurred atoms
 * ***************************************************************************
 */
void getMinMax(ATOM *atom_list, int size, float min[3], float max[3])
{
	int i, j;
	float maxRad=0.0;
	float tempRad;

	
	min[0] = max[0] = atom_list[0].x;
	min[1] = max[1] = atom_list[0].y;
	min[2] = max[2] = atom_list[0].z;
	maxRad = atom_list[0].radius * sqrt(1.+log(EPSILON)/BLOBBYNESS);
	
	for(j = 1; j < size; j++) {
		if (atom_list[j].x < min[0])
		  min[0] = atom_list[j].x;
		if (atom_list[j].y < min[1])
		  min[1] = atom_list[j].y;
		if (atom_list[j].z < min[2])
		  min[2] = atom_list[j].z;
		if (atom_list[j].x > max[0])
		  max[0] = atom_list[j].x;
		if (atom_list[j].y > max[1])
		  max[1] = atom_list[j].y;
		if (atom_list[j].z > max[2])
		  max[2] = atom_list[j].z;
		
		tempRad = atom_list[j].radius * sqrt(1.0 + log(EPSILON)/BLOBBYNESS);
		if (maxRad < tempRad)
		  maxRad = tempRad;
	}
	
	for(i = 0; i < 3; i++) {
	  min[i] -= maxRad;
	  max[i] += maxRad;
	}
}
