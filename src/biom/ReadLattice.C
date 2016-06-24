/*
 * ***************************************************************************
 * GAMER = < Geometry-preserving Adaptive MeshER >
 * Copyright (C) 2007-2010 -- Michael Holst and Zeyun Yu
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
 * File:     ReadLattice.C 
 *
 * Author:   Johan Hake (joha.hake@gmail.com)
 *
 * Purpose:  Provide functionality to read in Lattice data
 * ***************************************************************************
 */

// FIXME: This is not included by maloc?
#include <unistd.h>
#include <gamer/biom.h>

// Read a short integer, swapping the bytes
int read_short_int(FILE *fptr, short int* n)
{
   unsigned char *cptr, tmp;

   if (fread(n, sizeof(short int), 1, fptr) != 1)
      return(false);
   cptr = (unsigned char*)n;
   tmp = cptr[0];
   cptr[0] = cptr[1];
   cptr[1] =tmp;

   return(true);
}

// Read an unsigned integer from a char buffer, swapping the bytes
unsigned int get_uint(unsigned char*& cptr)
{
  unsigned int* n;
  char tmp;

  // Get the number
  n = (unsigned int*) cptr;

  // Swap bytes
  tmp = cptr[0];
  cptr[0] = cptr[3];
  cptr[3] = tmp;
  tmp = cptr[1];
  cptr[1] = cptr[2];
  cptr[2] = tmp;
  
  // Increase the char pointer
  cptr += 4;
  
  return *n;
}

// Read an integer from a char buffer, swapping the bytes
int get_int(unsigned char*& cptr)
{
  int* n;
  char tmp;

  // Get the number
  n = (int*) cptr;

  // Swap bytes
  tmp = cptr[0];
  cptr[0] = cptr[3];
  cptr[3] = tmp;
  tmp = cptr[1];
  cptr[1] = cptr[2];
  cptr[2] = tmp;
  
  // Increase the char pointer
  cptr += 4;
  
  return *n;
}

// Read a floating point number from a char buffer, swapping the bytes
// Assume IEEE format
float get_float(unsigned char*& cptr)
{
  float* n;
  char tmp;

  // Get the number
  n = (float*) cptr;

  // Swap bytes
  tmp = cptr[0];
  cptr[0] = cptr[3];
  cptr[3] = tmp;
  tmp = cptr[1];
  cptr[1] = cptr[2];
  cptr[2] = tmp;
  
  // Increase the char pointer
  cptr += 4;
  
  return *n;
}

// Read a double from a char buffer, swapping the bytes
// Assume IEEE
double get_double(unsigned char*& cptr)
{
  double* n;
  unsigned char tmp;

  // Get the number
  n = (double*) cptr;

  tmp = cptr[0];
  cptr[0] = cptr[7];
  cptr[7] = tmp;
  tmp = cptr[1];
  cptr[1] = cptr[6];
  cptr[6] = tmp;
  tmp = cptr[2];
  cptr[2] = cptr[5];
  cptr[5] =tmp;
  tmp = cptr[3];
  cptr[3] = cptr[4];
  cptr[4] = tmp;
  
  // Increase the char pointer
  cptr += 8;
  
  return *n;
} 

void load_lattice_file(char* input_name, unsigned int& xdim, 
		       unsigned int& ydim, unsigned int& zdim, float*& data, 
		       double*& dim_mat, bool padd, float padd_value)
{
  FILE* fptr;
  unsigned int file_size, i, j, k, read_version, write_version, smsize, dim;
  unsigned int data_type, compression_type, elem_size, dsize;
  size_t read_size, uint_size, double_size;
  unsigned char* buffer;
  unsigned char* cptr;
  char* header;

  // Open file
  if ((fptr=fopen(input_name, "rb"))==NULL) {fputs("File error\n", stderr); exit(ENOENT);}
  
  // Obtain file size:
  fseek (fptr, 0, SEEK_END);
  file_size = ftell(fptr);
  rewind (fptr);
  
  // Allocate memory to contain the whole file:
  buffer = (unsigned char*) malloc (sizeof(unsigned char)*file_size);
  if (buffer == NULL) {fputs("Memory error\n", stderr); exit(ENOMEM);}
  
  // Copy the file into the buffer and close file
  read_size = fread(buffer, 1, file_size, fptr);
  if (read_size != file_size) {fputs("Reading error\n", stderr); exit(EIO);}
  fclose(fptr);

  printf("File size: %d\n", file_size);

  header = (char*)buffer;

  // Check header
  if (strncmp(header, "LATTICE*IBT*KA", 14)){
    fputs("Not a Lattice file\n", stderr); 
    exit(EIO);
  };

  // Get the char pointer to the first data
  cptr = &buffer[14];

  // Read 'read_version'
  read_version = get_uint(cptr);
  printf("read version: %d\n", read_version);

  // Read 'write_version'
  write_version = get_uint(cptr);
  printf("write version: %d\n", write_version);

  // Read 'smsize'
  smsize = get_uint(cptr);
  printf("smsize: %d\n", smsize);

  // Read 'xdim'
  xdim = get_uint(cptr);
  printf("xdim: %d\n", xdim);

  // Read 'ydim'
  ydim = get_uint(cptr);
  printf("ydim: %d\n", ydim);

  // Read 'zdim'
  zdim = get_uint(cptr);
  printf("zdim: %d\n", zdim);

  // Read 'dim'
  dim = get_uint(cptr);
  printf("dim: %d\n", dim);

  // Read 'dim mat'
  dim_mat = (double*) malloc (sizeof(double)*4*dim);
  printf("Dimension matrix");
  for (i = 0; i < dim; i++){
    printf("\n");
    for (j = 0; j < 4; j++){
      dim_mat[i*4 + j] = get_double(cptr);
      printf("%f, ", dim_mat[i*4 + j]);
    }
  }
  
  printf("\n");
  // Read 'data_type'
  data_type = get_uint(cptr);
  printf("data_type: %d\n", data_type);

  // Read 'compression_type'
  compression_type = get_uint(cptr);
  printf("compression_type: %d\n", compression_type);

  if (compression_type){fputs("Compression type not supported\n", stderr); exit(EIO);}

  // Read 'elem_size'
  elem_size = get_uint(cptr);
  printf("elem_size: %d\n", elem_size);

  // Skip alignment
  cptr += 88; 

  // Read 'dsize'
  dsize = get_uint(cptr);
  printf("dsize: %d\n", dsize);

  // If padd create a padding around the data of size 1
  if (padd)
  {
    xdim+=2;
    ydim+=2;
    zdim+=2;
  }

  // Allocate memory
  data = (float*) malloc(sizeof(float)*xdim*ydim*zdim);

  // Read in data dependent on the data type
  switch (data_type) {
  case 0:
    
    if (padd)
    {
      // Initialize all data points with padd value
      for (i = 0; i < xdim*ydim*zdim; i++)
	data[i] = padd_value;

      for (k = 1; k < zdim-1; k++)
	for (j = 1; j < ydim-1; j++) 
	  for (i = 1; i < xdim-1; i++)
	    data[IndexVect(i,j,k)] = (*cptr++ == 255) ? 0.0 : 1.0;
    }
    else
      for (k = 0; k < zdim; k++)
	for (j = 0; j < ydim; j++) 
	  for (i = 0; i < xdim; i++)
	    data[IndexVect(i,j,k)] = (*cptr++ == 255) ? 0.0 : 1.0;
      
    break;
  case 4:

    if (padd)
    {
      // Initialize all data points with padd value
      for (i = 0; i < xdim*ydim*zdim; i++)
	data[i] = padd_value;
    
      // Fill inside of cell with data
      for (k = 1; k < zdim-1; k++)
	for (j = 1; j < ydim-1; j++) 
	  for (i = 1; i < xdim-1; i++)
	    data[IndexVect(i,j,k)] = get_float(cptr);
    }
    else
      for (k = 0; k < zdim; k++)
	for (j = 0; j < ydim; j++) 
	  for (i = 0; i < xdim; i++)
	    data[IndexVect(i,j,k)] = get_float(cptr);
    break;
  default:
    fputs("Data type not supported\n", stderr); exit(EIO);
  }

  // Release memory
  free(buffer);

}

/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_readLattice
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Read segmentation stored in a Lattice file and generate a 
 *           surface mesh
 * ***************************************************************************
 */
SurfaceMesh* SurfaceMesh_readLattice(char* segmentation_filename, 
				     float isovalue,
				     bool  padd)
{
  return SurfaceMesh_readLattice(segmentation_filename, NULL, 
				 isovalue, padd);
}


/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_readLattice
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Read segmentation stored in a Lattice file and generate a 
 *           surface mesh, potential using intensity values
 * ***************************************************************************
 */
SurfaceMesh* SurfaceMesh_readLattice(char* segmentation_filename, 
				     char* intensity_filename, 
				     float isovalue,
				     bool  padd)
{
  SPNT *holelist;
  unsigned int xdim0, ydim0, zdim0, xdim1, ydim1, zdim1;
  float* segmentation_data;
  float* intensity_data;
  double* mat_dim = NULL;
  double* mat_dim_intensity = NULL;
  unsigned int basename_index;
  char* basename_end_pointer;
  SurfaceMesh* surfmesh;

  printf("Load Lattice\n");

  // Check if input is a .lat file
  basename_end_pointer = strstr(segmentation_filename, ".lat");
  if (basename_end_pointer == NULL){
    printf("No \'.lat\' file\n");
    return NULL;
  }

  // Check if file exists and are readable
  if (access(segmentation_filename, R_OK)){
    printf("Not readable\n");
    return NULL;
  }
  
  // Load the segmented data
  load_lattice_file(segmentation_filename, xdim0, ydim0, zdim0, 
		    segmentation_data, mat_dim, padd, 0.0);
  
  // Check if we intensity_filename is given
  if (intensity_filename != NULL)
  {

    printf("Intensity file given\n");

    // Check if input is a .flat file
    basename_end_pointer = strstr(intensity_filename, ".flat");
    if (basename_end_pointer == NULL){
      printf("No \'.flat\' file\n");
      return NULL;
    }
    
    // Check if file exists and are readable
    if (access(intensity_filename, R_OK)){
      printf("Not readable\n");
      return NULL;
    }
  
    // Load the intensity data
    load_lattice_file(intensity_filename, xdim1, ydim1, zdim1, intensity_data, 
		      mat_dim_intensity, padd, -1.);
   
    // Check if dimensions fit
    if (xdim0 != xdim1 && ydim0 != ydim1 && zdim0 != zdim1)
    {
      printf("Dimension does not fit\n");
      return NULL;
    }

    // Call the marhching cube rutine with intensity data
    surfmesh = SurfaceMesh_marchingCube(xdim0, ydim0, zdim0, segmentation_data, 0.5,
					intensity_data, isovalue, &holelist);
    
    // Centralize and scale data
    SurfaceMesh_centeralize(surfmesh);
    SurfaceMesh_scale(surfmesh, mat_dim[0], mat_dim[5], mat_dim[10]);

    // Free memory
    free(segmentation_data);
    free(intensity_data);
    free(mat_dim);
    free(mat_dim_intensity);
  
  }
  else
  {

    // Call the marhching cube rutine without intensity data
    surfmesh = SurfaceMesh_marchingCube(xdim0, ydim0, zdim0, segmentation_data, 
					isovalue, &holelist);
  
    // Centralize and scale data
    SurfaceMesh_centeralize(surfmesh);
    SurfaceMesh_scale(surfmesh, mat_dim[0], mat_dim[5], mat_dim[10]);

    // Free memory
    free(segmentation_data);
    free(mat_dim);

  }

  // Return surface mesh
  return surfmesh;
}

