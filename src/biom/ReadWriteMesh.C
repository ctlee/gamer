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
 * Routine:  Sort cell connectivities (Needed by DOLFIN)
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  swap the bytes when LITTLE_ENDIAN is enabled
 * ***************************************************************************
 */
#include <vector>
#include <limits>

void sort_connectivity(const std::vector<unsigned int>& conn, unsigned int* indices)
{
  unsigned int min_val=std::numeric_limits<unsigned int>::max(), max_val=0;
  unsigned int min_ind=0, max_ind=3;
  bool left_overs[4] = {true, true, true, true};

  // First find min and max indices
  for (unsigned int i=0; i<4; i++)
  {
    if (conn[i]<min_val)
    {
      min_val = conn[i];
      min_ind = i;
    }
    
    if (conn[i]>max_val)
    {
      max_val = conn[i];
      max_ind = i;
    }
  }
  
  // Store min and max indices
  indices[0] = min_ind;
  indices[3] = max_ind;

  // Sort the two other indices
  left_overs[min_ind] = false;
  left_overs[max_ind] = false;
  
  // Find the other indices
  unsigned int other_inds[2] = {0, 0};
  unsigned int ind_into_other_inds = 0;

  for (unsigned int i=0; i<4; i++)
  {
    if (left_overs[i])
      other_inds[ind_into_other_inds++] = i;
    if (ind_into_other_inds>1)
      break;
  } 

  if (conn[other_inds[0]]<conn[other_inds[1]])
  {
    indices[1] = other_inds[0];
    indices[2] = other_inds[1];
  }
  else
  {
    indices[1] = other_inds[1];
    indices[2] = other_inds[0];
  }
}

/* 
 * ***************************************************************************
 * File:     ReadWriteMesh.C
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Some input and output routines
 * ***************************************************************************
 */

#include <stdexcept>
#include <vector>
#include <set>
#include <map>
#include <string>

#include <gamer/biom.h>
#include <gamer/tetgen.h>
#include "gamercf.h"

/* 
 * All ISO-C headers are coming to us via MALOC, through the includes above.
 * We now include all architecture-dependent headers.
 * This buries the architecture-dependent stuff in library, not headers.
 */

#if defined(HAVE_UNISTD_H)
#   include <unistd.h>
#endif

#if defined(HAVE_SYS_TYPES_H)
#   include <sys/types.h>
#endif

#if defined(HAVE_SYS_TIME_H)
#   include <sys/time.h>
#endif

#if defined(HAVE_SYS_TIMES_H)
#   include <sys/times.h>
#endif

#if defined(HAVE_SYS_STAT_H)
#   include <sys/stat.h>
#endif


/*
 * ***************************************************************************
 * Routine:  swap_buffer
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  swap the bytes when LITTLE_ENDIAN is enabled
 * ***************************************************************************
 */
void swap_buffer(char *buffer, int count, int typesize)
{
  char sbuf[4];
  int i;
  int temp = 1;
  unsigned char* chartempf = (unsigned char*) &temp;
  if(chartempf[0] > '\0') {
  
    // swapping isn't necessary on single byte data
    if (typesize == 1)
      return;
  
  
    for (i=0; i < count; i++)
      {
	memcpy(sbuf, buffer+(i*typesize), typesize);
      
	switch (typesize)
	  {
	  case 2:
	    {
	      buffer[i*typesize] = sbuf[1];
	      buffer[i*typesize+1] = sbuf[0];
	      break;
	    }
	  case 4:
	    {
	      buffer[i*typesize] = sbuf[3];
	      buffer[i*typesize+1] = sbuf[2];
	      buffer[i*typesize+2] = sbuf[1];
	      buffer[i*typesize+3] = sbuf[0];
	      break;
	    }
	  default:
	    break;
	  }
      }

  }
}



/*
 * ***************************************************************************
 * Routine:  ReadRawiv
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Read the volume in rawiv format 
             (see Chandrajit Bajaj's group for details of this format)
 * ***************************************************************************
 */
void ReadRawiv(int *xd, int *yd, int *zd, float **dataset, char *input_name,
               float *span_t, float *orig_t)
{
  float c_float;
  unsigned char c_unchar;
  unsigned short c_unshort;
  int i,j,k;
  float *data;
  int xdim,ydim,zdim;
  float maxraw;
  float minraw;
  float minext[3], maxext[3];
  int nverts, ncells;
  unsigned int dim[3];
  float orig[3], span[3];

  struct stat filestat;
  size_t size[3];
  int datatype;
  int found;
  FILE *fp;
  size_t ret;
  

  if ((fp=fopen(input_name, "rb"))==NULL){
    printf("read error...\n");
    exit(0);
  }
  stat(input_name, &filestat);

  /* reading RAWIV header */
  ret = fread(minext, sizeof(float), 3, fp);
  ret = fread(maxext, sizeof(float), 3, fp);
  ret = fread(&nverts, sizeof(int), 1, fp);
  ret = fread(&ncells, sizeof(int), 1, fp);
#ifdef _LITTLE_ENDIAN
  swap_buffer((char *)minext, 3, sizeof(float));
  swap_buffer((char *)maxext, 3, sizeof(float));
  swap_buffer((char *)&nverts, 1, sizeof(int));
  swap_buffer((char *)&ncells, 1, sizeof(int));
#endif  

  size[0] = 12 * sizeof(float) + 2 * sizeof(int) + 3 * sizeof(unsigned int) +
    nverts * sizeof(unsigned char);
  size[1] = 12 * sizeof(float) + 2 * sizeof(int) + 3 * sizeof(unsigned int) +
    nverts * sizeof(unsigned short);
  size[2] = 12 * sizeof(float) + 2 * sizeof(int) + 3 * sizeof(unsigned int) +
    nverts * sizeof(float);
  
  found = 0;
  for (i = 0; i < 3; i++)
    if (size[i] == (unsigned int)filestat.st_size)
      {
        if (found == 0)
          {
            datatype = i;
            found = 1;
          }
      }
  if (found == 0)
    {
      printf("Corrupted file or unsupported dataset type\n");
      exit(5);
    }
 
  ret = fread(dim, sizeof(unsigned int), 3, fp);
  ret = fread(orig, sizeof(float), 3, fp);
  ret = fread(span, sizeof(float), 3, fp);
#ifdef _LITTLE_ENDIAN
  swap_buffer((char *)dim, 3, sizeof(unsigned int));
  swap_buffer((char *)orig, 3, sizeof(float));
  swap_buffer((char *)span, 3, sizeof(float));
#endif 

  span_t[0] = span[0];
  span_t[1] = span[1];
  span_t[2] = span[2];
  orig_t[0] = orig[0];
  orig_t[1] = orig[1];
  orig_t[2] = orig[2];

  xdim = dim[0];
  ydim = dim[1];
  zdim = dim[2];
  
  data = (float *)malloc(sizeof(float)*xdim*ydim*zdim);
  
  maxraw = -999999.f;
  minraw = 999999.f;
  if (datatype == 0) {
    printf("data type: unsigned char \n");
    for (i=0; i<zdim; i++)
      for (j=0; j<ydim; j++)
        for (k=0; k<xdim; k++) {
          ret = fread(&c_unchar, sizeof(unsigned char), 1, fp);
          data[IndexVect(k,j,i)]=(float)c_unchar;
	  if (c_unchar > maxraw)
            maxraw = c_unchar;
          if (c_unchar < minraw)
            minraw = c_unchar;
        }
  }
  else if (datatype == 1) {
    printf("data type: unsigned short \n");
    for (i=0; i<zdim; i++)
      for (j=0; j<ydim; j++)
	for (k=0; k<xdim; k++) {
	  ret = fread(&c_unshort, sizeof(unsigned short), 1, fp);
#ifdef _LITTLE_ENDIAN
	  swap_buffer((char *)&c_unshort, 1, sizeof(unsigned short));
#endif 
	  data[IndexVect(k,j,i)]=(float)c_unshort;
	  
	  if (c_unshort > maxraw)
	    maxraw = c_unshort;
	  if (c_unshort < minraw)
	    minraw = c_unshort;
	}
  }
  else if (datatype == 2) {
    printf("data type: float \n");
    for (i=0; i<zdim; i++) 
      for (j=0; j<ydim; j++)
	for (k=0; k<xdim; k++) {
	  ret = fread(&c_float, sizeof(float), 1, fp);
#ifdef _LITTLE_ENDIAN
	  swap_buffer((char *)&c_float, 1, sizeof(float));
#endif 
	  data[IndexVect(k,j,i)]=c_float;
            
	  if (c_float > maxraw)
	    maxraw = c_float;
	  if (c_float < minraw)
	    minraw = c_float;
	}
  }
  else {
    printf("error\n");
    fclose(fp);
    exit(1);
  }
  
  fclose(fp);
  printf("minimum = %f,   maximum = %f \n",minraw,maxraw);

  for (i=0; i<zdim; i++) 
    for (j=0; j<ydim; j++)
      for (k=0; k<xdim; k++) {
        data[IndexVect(k,j,i)] = 255*(data[IndexVect(k,j,i)] - 
                                         minraw)/(maxraw-minraw); 
      }
  
  printf("dimension: %d X %d X %d\n",xdim,ydim,zdim);
  *xd = xdim;
  *yd = ydim;
  *zd = zdim;
  *dataset = data;

}



/*
 * ***************************************************************************
 * SubRoutine:  rgb_to_marker
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Transform an rgb value to an integer
 *           
 * ***************************************************************************
 */
int rgb_to_marker(float r, float g, float b)
{
  if (r < 0 | r > 1 | g < 0 | g > 1 | b < 0 | b > 1 ){
    printf("Expected individual RGB value to be betwen 0 and 1.\n");
    exit(1);
  }
  
  return round(r*10)*121 + round(g*10)*11 + round(b*10);
}


/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_readPoly
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Read a tetgen poly file representing a Piecewise Linear Complex
 * ***************************************************************************
 */
SurfaceMesh* SurfaceMesh_readPoly(char *input_name)
{
  int n, i, j;
  
  // Use Tetgen's tetgenio to load poly file
  tetgenio tetio;
  
  if (!tetio.load_poly(input_name))
  {
    printf("Read error. File \'%s\' could not be read.\n", input_name);
    return NULL;
  }
  
  // Create the SurfaceMesh
  SurfaceMesh* surfmesh = SurfaceMesh_ctor(tetio.numberofpoints, tetio.numberoffacets);
  
  // Get vertex information
  for (n=0; n<tetio.numberofpoints; n++)
  {
    surfmesh->vertex[n].x = tetio.pointlist[n*3];
    surfmesh->vertex[n].y = tetio.pointlist[n*3+1];
    surfmesh->vertex[n].z = tetio.pointlist[n*3+2];
    surfmesh->vertex[n].sel = true;
    if (tetio.pointmarkerlist!=NULL)
      surfmesh->vertex[n].m = tetio.pointmarkerlist[n];
    
  }
  
  // Get all face informations
  for (n=0; n<tetio.numberoffacets; n++)
  {

    // Check that we have something sane
    if (tetio.facetlist[n].numberofpolygons != 1)
    {
      printf("Read error: %s. A Facet can only include one polygon.\n", input_name);
      return NULL;
    }
    
    if (tetio.facetlist[n].polygonlist[0].numberofvertices!=3)
    {
      printf("Read error: %s. Facets can only be triangles.\n", input_name);
      return NULL;
    }
    
    // Assign vertex indices
    surfmesh->face[n].a = tetio.facetlist[n].polygonlist[0].vertexlist[0]-1;
    surfmesh->face[n].b = tetio.facetlist[n].polygonlist[0].vertexlist[1]-1;
    surfmesh->face[n].c = tetio.facetlist[n].polygonlist[0].vertexlist[2]-1;
    surfmesh->face[n].sel = true;

    // Check for facet marker info
    if (tetio.facetmarkerlist!=NULL)
      surfmesh->face[n].m = tetio.facetmarkerlist[n];
    
  }
  
  // Return the mesh
  return surfmesh;

}

/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_writePoly
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Write a tetgen poly file representing a Piecewise Linear Complex
 * ***************************************************************************
 */
void SurfaceMesh_writePoly(SurfaceMesh* surfmesh, char *filename)
{
  int i;
  FILE* fout;

  // Try to open the file
  if ((fout=fopen(filename, "wb"))==NULL){
    printf("write error...\n");
    return;
  };    

  fprintf(fout, "# Part 1 - node list\n");
  fprintf(fout, "%d 3 0 0\n", surfmesh->num_vertices);
  for (i = 0; i < surfmesh->num_vertices; i++) 
    fprintf(fout, "%d %17.10e    %17.10e    %17.10e\n", i+1, surfmesh->vertex[i].x, \
	    surfmesh->vertex[i].y, surfmesh->vertex[i].z);
  
  fprintf(fout, "# Part 2 - facet list\n");
  fprintf(fout, "%d 1\n", surfmesh->num_vertices);
  for (i = 0; i < surfmesh->num_faces; i++) 
  {
    fprintf(fout, "1 0 %d\n", surfmesh->face[i].m);
    fprintf(fout, "3 %d %d %d\n",surfmesh->face[i].a+1, surfmesh->face[i].b+1, \
	    surfmesh->face[i].c+1);
  }

  fprintf(fout, "# Part 3 - hole list\n0\n");

  fclose(fout);
}

/*
 * ***************************************************************************
 * Routine:  SurfaceMesh_readOFF
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Read a surface mesh in OFF format, if a volumetric mesh is 
 *           provided, a surface mesh will be extracted
 * ***************************************************************************
 */
SurfaceMesh* SurfaceMesh_readOFF(char *input_name)
{
  unsigned int n,m;
  unsigned int a,b,c,d;
  float x,y,z, color_r, color_g, color_b, color_a;
  TeTraMesh *volmesh = NULL;
  int v_n, t_n, e_n, character;
  char line[256];
  FILE *fin = NULL;
  fpos_t position;
  bool has_marker = false;
  SurfaceMesh* surfmesh = NULL;
  

  if ((fin = fopen(input_name, "r"))==NULL){
    printf("Read error. File \'%s\' could not be read.\n", input_name);
    return NULL;
  };
  
  while (fgets(line, 256, fin) != NULL) {
    if (line[0]=='O' && line[1]=='F' && line[2]=='F')
      break;
  }
  
  if (fscanf(fin,"%d %d %d\n", &v_n, &t_n, &e_n) != 3){
    printf("Read error. Expected 3 integer number for the number of vertices "
	   "and simplices.\n");
    return NULL;
  }
    
  printf("   vertices: %d --- simplices: %d \n", v_n, t_n);
  
  // Assign memory
  surfmesh = SurfaceMesh_ctor(v_n, t_n);
  
  // Read vertex coordinates
  for (n = 0; n < surfmesh->num_vertices; n++) {
    if (fscanf(fin, "%f %f %f\n", &x, &y, &z) != 3){
      printf("Read error. Expected 3 floats for the coordinates of vertex: %d.\n", n);
      return NULL;
    }

    surfmesh->vertex[n].x = x;
    surfmesh->vertex[n].y = y;
    surfmesh->vertex[n].z = z;
  }
  
  // Read the number of vertices per simplex from the first line of simplices
  n = fscanf(fin, "%d", &m);
    
  // Check format of the read simplices
  if (m != 3 && m != 4) {
    printf("Read error. Expected a 3 or 4 for the first value in the first simplex line.\n");
    return NULL;
  }
  
  // Input is surface mesh
  else if (m == 3) {
    printf("   Input is surface mesh.\n");
    
    // Read the rest of the first line
    if (fscanf(fin, "%d %d %d", &a, &b, &c) != 3){
      printf("Read error. Expected 3 integers for the first simplex.\n");
      return NULL;
    }
    
    // Get position of file and check if we have more symbols
    fgetpos(fin, &position);

    // Read any blanks
    character = fgetc(fin);
    while (character == ' '){
      character = fgetc(fin);
    }

    // If we do not have en end of line character we have marker
    if (character != '\n'){

      // Set marker flag
      has_marker = true;

      // Rewind the file pointer
      fsetpos(fin, &position);

      // Read first rgb values
      if (fscanf(fin, "%f %f %f %f", &color_r, &color_g, &color_b, &color_a) != 4){
	printf("Read error. Expected 4 floats for the RGBA values of the first face.\n");
	return NULL;
      }
      
      // Get the first face marker
      surfmesh->face[0].m = rgb_to_marker(color_r, color_g, color_b);

      while (fgetc(fin) != '\n'){}
    }

    // Assign first simplex
    surfmesh->face[0].a = a;
    surfmesh->face[0].b = b;
    surfmesh->face[0].c = c;
    
    // Read the rest of the simplices
    for (n = 1; n < surfmesh->num_faces; n++) {
      if (fscanf(fin, "%d %d %d %d", &m, &a, &b, &c) != 4){
	printf("Read error. Expected 4 integers for simplex %d.\n", n);
	return NULL;
      }
      surfmesh->face[n].a = a;
      surfmesh->face[n].b = b;
      surfmesh->face[n].c = c;

      // If we have face markers
      if (has_marker){
      
	// Read first rgb values
	if (fscanf(fin, "%f %f %f %f", &color_r, &color_g, &color_b, &color_a) != 4){
	  printf("Read error. Expected 4 floats for the RGBA values of face: %d.\n", n);
	  return NULL;
	}
      
	// Get the other face markers
	surfmesh->face[n].m = rgb_to_marker(color_r, color_g, color_b);

      }

      // Skip any additional character on this line
      while (fgetc(fin) != '\n'){}
    }

    // Close file
    fclose(fin);
  }

  // Input is tetrahedral mesh.
  else {
    printf("   Input is volumetric mesh...\n");
    fclose(fin);

    // Create tetrahedal mesh
    volmesh = (TeTraMesh*)malloc(sizeof(TeTraMesh));
    volmesh->num_vertices = surfmesh->num_vertices;
    volmesh->num_cells = surfmesh->num_faces;
    volmesh->neighbor = NULL;

    // We have allready read the vertex coordinates
    volmesh->vertex = surfmesh->vertex;

    // Assign memory for faces
    volmesh->face = (INT4VECT *)malloc(sizeof(INT4VECT)*volmesh->num_cells);

    // Read the rest of the first line
    if (fscanf(fin, "%d %d %d %d", &a, &b, &c, &d) != 4)
    {
      printf("Read error. Expected 4 integers for the first simplex.\n");
      return NULL;
    }
          
    // Skip any additional character on this line
    while (fgetc(fin) != '\n'){}

    // Assign first simplex from allready read data
    volmesh->face[0].a = a;
    volmesh->face[0].b = b;
    volmesh->face[0].c = c;
    volmesh->face[0].d = d;
    
    // Read the rest of the simplices
    for (n = 1; n < volmesh->num_cells; n++)
    {
      if (fscanf(fin, "%d %d %d %d %d", &m, &a, &b, &c, &d) != 5)
      {
	printf("Read error. Expected 5 integers for simplex %d.\n", n);
	return NULL;
      }
      volmesh->face[n].a = a;
      volmesh->face[n].b = b;
      volmesh->face[n].c = c;
      volmesh->face[n].d = d;

      // Skip any additional character on this line
      while (fgetc(fin) != '\n'){}
    }
    
    fclose(fin);

    SurfaceExtract(volmesh, surfmesh);
    free(volmesh->vertex);
    free(volmesh->face);
    free(volmesh);
    
  }

  return surfmesh;
}




/*
 * ***************************************************************************
 * Routine:  write_rawiv_float
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  write the volume into a rawiv file
 * ***************************************************************************
 */
void write_rawiv_float(FILE* fp, float *dataset, int dim[3], 
		       float min[3], float max[3])
{
  int i, j, k;
  float c_float;
  float orig[3],span[3];
  int xdim,ydim,zdim;
  int nverts,ncells;
  
  
  xdim = dim[0];
  ydim = dim[1];
  zdim = dim[2];
    
  orig[0] = min[0];
  orig[1] = min[1];
  orig[2] = min[2];
  span[0] = (max[0] - min[0]) / (float)(xdim-1);
  span[1] = (max[1] - min[1]) / (float)(ydim-1);
  span[2] = (max[2] - min[2]) / (float)(zdim-1);
  nverts = xdim*ydim*zdim;
  ncells = (xdim-1)*(ydim-1)*(zdim-1);
  

#ifdef _LITTLE_ENDIAN
  swap_buffer((char *)min, 3, sizeof(float));
  swap_buffer((char *)max, 3, sizeof(float));
  swap_buffer((char *)&nverts, 1, sizeof(int));
  swap_buffer((char *)&ncells, 1, sizeof(int));
  swap_buffer((char *)dim, 3, sizeof(unsigned int));
  swap_buffer((char *)orig, 3, sizeof(float));
  swap_buffer((char *)span, 3, sizeof(float));
#endif 
  
  fwrite(min, sizeof(float), 3, fp);
  fwrite(max, sizeof(float), 3, fp);
  fwrite(&nverts, sizeof(int), 1, fp);
  fwrite(&ncells, sizeof(int), 1, fp);
  fwrite(dim, sizeof(unsigned int), 3, fp);
  fwrite(orig, sizeof(float), 3, fp);
  fwrite(span, sizeof(float), 3, fp);
#ifdef _LITTLE_ENDIAN
  swap_buffer((char *)min, 3, sizeof(float));
  swap_buffer((char *)max, 3, sizeof(float));
  swap_buffer((char *)&nverts, 1, sizeof(int));
  swap_buffer((char *)&ncells, 1, sizeof(int));
  swap_buffer((char *)dim, 3, sizeof(unsigned int));
  swap_buffer((char *)orig, 3, sizeof(float));
  swap_buffer((char *)span, 3, sizeof(float));
#endif 

  for (k=0; k<zdim; k++) 
    for (j=0; j<ydim; j++)
      for (i=0; i<xdim; i++) {
        c_float = dataset[IndexVect(i,j,k)];
#ifdef _LITTLE_ENDIAN
        swap_buffer((char *)&c_float, 1, sizeof(float));
#endif 
        fwrite(&c_float, sizeof(float), 1, fp);
#ifdef _LITTLE_ENDIAN
        swap_buffer((char *)&c_float, 1, sizeof(float));
#endif 
      }
  
  fclose(fp);
}


/*
 * ***************************************************************************
 * Routine:  WriteSurfaceMeshOFF
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  write a surface mesh into an OFF file
 * ***************************************************************************
 */
void SurfaceMesh_writeOFF(SurfaceMesh *surfmesh, char* filename)
{
  int i;
  FILE* fout;

  // Try to open the file
  if ((fout=fopen(filename, "wb"))==NULL){
    printf("write error...\n");
    return; 
  };    

  fprintf(fout, "OFF\n");
  fprintf(fout, "%d %d %d\n",surfmesh->num_vertices,surfmesh->num_faces,surfmesh->num_vertices+surfmesh->num_faces-2);

  for (i = 0; i < surfmesh->num_vertices; i++) 
    fprintf(fout, "%17.10e    %17.10e    %17.10e\n",surfmesh->vertex[i].x,surfmesh->vertex[i].y,surfmesh->vertex[i].z);
  
  for (i = 0; i < surfmesh->num_faces; i++) {
    fprintf(fout, "3 %d %d %d\n",surfmesh->face[i].a,surfmesh->face[i].b,surfmesh->face[i].c);
  }
  fclose(fout);
}

/*
 * ***************************************************************************
 * Routine:  write_off
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Write a surface mesh into an OFF format
 * ***************************************************************************
 */
void GemMesh_writeOFF(GemMesh* gem_mesh, char* filename)
{
  long i;
  
  FILE* fout;

  // Try to open the file
  if ((fout=fopen(filename, "wb"))==NULL){
    printf("write error...\n");
    return; 
  };    

  fprintf(fout, "OFF\n");
  fprintf(fout, "%d  %d  0\n", gem_mesh->num_vertices, gem_mesh->num_cells);
  for (i = 0; i < gem_mesh->num_vertices; i++) 
    fprintf(fout, "%.12f %.12f %.12f\n", gem_mesh->vv[i].x, gem_mesh->vv[i].y, gem_mesh->vv[i].z);
    
  for (i = 0; i < gem_mesh->num_cells; i++) 
    fprintf(fout, "4  %8d  %8d  %8d  %8d\n", gem_mesh->ss[i].na, gem_mesh->ss[i].nb, 
	    gem_mesh->ss[i].nc, gem_mesh->ss[i].nd);

  fclose(fout);
}

/*
 * ***************************************************************************
 * Routine:  GemMesh_writeDolfin
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Write tetrahedralized mesh to disk (DOLFIN formats)
 * ***************************************************************************
 */

void GemMesh_writeDolfin(GemMesh* out, char* filename)
{

  if (out->higher_order)
  {
    printf("Dolfin output does not support higher order mesh...\n");
    return; 
  }

  unsigned int i, j;
  FILE *fout;
  FETK_VX* vv;
  FETK_SS* ss;

  // Open file
  if ((fout=fopen(filename, "wb"))==NULL){
    printf("write error...\n");
    return; 
  };    

  fprintf(fout, "<?xml version=\"1.0\"?>\n");
  fprintf(fout, "<dolfin xmlns:dolfin=\"http://fenicsproject.org\">\n");
  fprintf(fout, "  <mesh celltype=\"tetrahedron\" dim=\"3\">\n");
  fprintf(fout, "    <vertices size=\"%d\">\n", out->num_vertices);

  // Write all vertices to file
  for (i = 0; i < out->num_vertices; i++)
  {
    vv = out->vv + i;
    fprintf(fout, "      <vertex index=\"%d\" x=\"%f\" y=\"%f\" z=\"%f\" />\n", i, 
	    vv->x, vv->y, vv->z);
  }

  fprintf(fout, "    </vertices>\n");
  fprintf(fout, "    <cells size=\"%d\">\n", out->num_cells);

  // Write all cells to file 
  unsigned int num_face_markers=0;
  std::vector<unsigned int> conn(4);
  std::vector<unsigned int> indices(4*out->num_cells);
  for (i = 0; i < out->num_cells; i++) 
  {
    
    ss = out->ss + i;

    // Sort vertices
    conn[0]=ss->na; conn[1]=ss->nb; conn[2]=ss->nc; conn[3]=ss->nd;
    sort_connectivity(conn, &indices[i*4]);

    fprintf(fout, "      <tetrahedron index=\"%d\" v0=\"%d\" v1=\"%d\" v2=\"%d\" "\
	    "v3=\"%d\" />\n", i, conn[indices[i*4]], conn[indices[i*4+1]], \
	    conn[indices[i*4+2]], conn[indices[i*4+3]]);

    // Count face markers
    if (ss->fa!=0)
      num_face_markers += 1;
    if (ss->fb!=0)
      num_face_markers += 1;
    if (ss->fc!=0)
      num_face_markers += 1;
    if (ss->fd!=0)
      num_face_markers += 1;
  }

  fprintf(fout, "    </cells>\n");
  fprintf(fout, "    <domains>\n");
  
  fprintf(fout, "      <mesh_value_collection name=\"m\" type=\"uint\" "
	  "dim=\"2\" size=\"%d\">\n", num_face_markers);

  // Iterate over cells and collect face markers which are different from 0
  for (i = 0; i < out->num_cells; i++) 
  {
    ss = out->ss + i;
    conn[0]=ss->fa; conn[1]=ss->fb; conn[2]=ss->fc; conn[3]=ss->fd;
    for (j=0; j<4; j++)
    {
      if (conn[indices[i*4+j]] != 0)
	fprintf(fout, "        <value cell_index=\"%d\" local_entity=\"%d\" "
		"value=\"%d\" />\n", i, j, conn[indices[i*4+j]]);
    }
  }

  fprintf(fout, "      </mesh_value_collection>\n");
  fprintf(fout, "      <mesh_value_collection name=\"m\" type=\"uint\" "
	  "dim=\"3\" size=\"%d\">\n", out->num_cells);

  // Iterate over cells and write cell markers
  for (i = 0; i < out->num_cells; i++) 
  {
    ss = out->ss + i;
    fprintf(fout, "        <value cell_index=\"%d\" local_entity=\"0\" "
	    "value=\"%d\" />\n", i, ss->mat);
  }

  fprintf(fout, "      </mesh_value_collection>\n");
  fprintf(fout, "    </domains>\n");
  fprintf(fout, "  </mesh>\n");
  fprintf(fout, "</dolfin>\n");
  
  fclose(fout);

}

/*
 * ***************************************************************************
 * Routine:  GemMesh_writeCarp
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Write tetrahedralized mesh to disk (Carp format)
 * ***************************************************************************
 */
void GemMesh_writeCarp(GemMesh* out, char* basename, int num_boundaries, 
		       char* boundaries[])
{

  if (out->higher_order)
  {
    printf("Carp output does not support higher order mesh...\n");
    return; 
  }

  char filename [256];
  unsigned int i;
  FILE *fout;
  FETK_VX* vv;
  FETK_SS* ss;
  INT3VECT* triangle;

  // Collect face information
  std::map<int, std::set<int> > vertex_markers;
  std::set<int> all_markers;
  
  for (i = 0; i < out->num_cells; i++)
  {
    ss = out->ss + i;

    // Count face markers
    if (ss->fa!=0)
    {
      if (all_markers.count(ss->fa)==0)
      {
	std::set<int> _vertices;
	all_markers.insert(ss->fa);
	vertex_markers[ss->fa] = _vertices;
      }
      vertex_markers[ss->fa].insert(ss->nb);
      vertex_markers[ss->fa].insert(ss->nc);
      vertex_markers[ss->fa].insert(ss->nd);
    }
    if (ss->fb!=0)
    {
      if (all_markers.count(ss->fb)==0)
      {
	std::set<int> _vertices;
	all_markers.insert(ss->fb);
	vertex_markers[ss->fb] = _vertices;
      }
      vertex_markers[ss->fb].insert(ss->na);
      vertex_markers[ss->fb].insert(ss->nc);
      vertex_markers[ss->fb].insert(ss->nd);
    }
    if (ss->fc!=0)
    {
      if (all_markers.count(ss->fc)==0)
      {
	std::set<int> _vertices;
	all_markers.insert(ss->fc);
	vertex_markers[ss->fc] = _vertices;
      }
      vertex_markers[ss->fc].insert(ss->na);
      vertex_markers[ss->fc].insert(ss->nb);
      vertex_markers[ss->fc].insert(ss->nd);
    }
    if (ss->fd!=0)
    {
      if (all_markers.count(ss->fd)==0)
      {
	std::set<int> _vertices;
	all_markers.insert(ss->fd);
	vertex_markers[ss->fd] = _vertices;
      }
      vertex_markers[ss->fd].insert(ss->na);
      vertex_markers[ss->fd].insert(ss->nb);
      vertex_markers[ss->fd].insert(ss->nc);
    }
  }

  // Check that number of facet markers equals the number passed boundary names
  if (all_markers.size()!=num_boundaries)
  {
    throw std::runtime_error("Number of passed boundaries does not match number "
			    "of facet markers in mesh."); 
  }

  // Output vertex file
  sprintf (filename, "%s.pts", basename);
  if ((fout=fopen(filename, "wb"))==NULL){
    printf("write error...\n");
    return; 
  };    

  fprintf(fout, "%d\n", out->num_vertices);
  for (i = 0; i < out->num_vertices; i++)
  {
    vv = out->vv + i;
    fprintf(fout, "%e, %e, %e\n", vv->x, vv->y, vv->z);
  }
  fclose(fout);

  // Output tet file
  sprintf (filename, "%s.elem", basename);
  if ((fout=fopen(filename, "wb"))==NULL){
    printf("write error...\n");
    return; 
  };    

  fprintf(fout, "%d\n", out->num_cells);
  for (i = 0; i < out->num_cells; i++)
  {
    ss = out->ss + i;

    // Check sign of determinant det(a - b, b - c, c - d)|
    // If negative change or of last two vertices
    //const double ax = out->vv[ss->na].x;
    //const double ay = out->vv[ss->na].y;
    //const double az = out->vv[ss->na].z;
    //
    //const double bx = out->vv[ss->nb].x;
    //const double by = out->vv[ss->nb].y;
    //const double bz = out->vv[ss->nb].z;
    //
    //const double cx = out->vv[ss->nc].x;
    //const double cy = out->vv[ss->nc].y;
    //const double cz = out->vv[ss->nc].z;
    //
    //const double dx = out->vv[ss->nd].x;
    //const double dy = out->vv[ss->nd].y;
    //const double dz = out->vv[ss->nd].z;

    /*
    // Formula for the determinant of tetrahedron from http://mathworld.wolfram.com
    const double det_j = (ax*(by*cz + dy*bz + cy*dz - cy*bz - by*dz - dy*cz) -
		          bx*(ay*cz + dy*az + cy*dz - cy*az - ay*dz - dy*cz) +
		          cx*(ay*bz + dy*az + by*dz - by*az - ay*dz - dy*bz) -
		          dx*(ay*bz + by*cz + cy*az - by*az - cy*bz - ay*cz));
    */
 
    //const double j11 = bx - ax; 
    //const double j12 = by - ay; 
    //const double j13 = bz - az;
    //const double j21 = cx - ax; 
    //const double j22 = cy - ay; 
    //const double j23 = cz - az;
    //const double j31 = dx - ax; 
    //const double j32 = dy - ay; 
    //const double j33 = dz - az;
    //
    //const double a11 = j22*j33 - j23*j32;  
    //const double a12 = j13*j32 - j12*j33;  
    //const double a13 = j12*j23 - j13*j22;
    //const double a21 = j23*j31 - j21*j33;  
    //const double a22 = j11*j33 - j13*j31;  
    //const double a23 = j13*j21 - j11*j23;
    //const double a31 = j21*j32 - j22*j31;  
    //const double a32 = j12*j31 - j11*j32;  
    //const double a33 = j11*j22 - j12*j21;
    //
    //const double det_j = j11*a11 + j12*a21 + j13*a31;

    // NOTE: Diffpack use 1 as its first cell number, also note that each 
    // NOTE: vertex index also need to be bumped with 1
    //if (det_j > 0)
    //  fprintf(fout, "Tt %d %d %d %d %d\n",		\
    //	      ss->na, ss->nb, ss->nc, ss->nd, ss->mat);
    //else
    //  fprintf(fout, "Tt %d %d %d %d %d\n",		\
    //	      ss->na, ss->nb, ss->nd, ss->nc, ss->mat);
    fprintf(fout, "Tt %d %d %d %d %d\n",    
	    ss->na, ss->nb, ss->nc, ss->nd, ss->mat);
  }

  fclose(fout);
  
  // Write surf file (the boundary)
  sprintf (filename, "%s.surf", basename);
  if ((fout=fopen(filename, "wb"))==NULL){
    printf("write error...\n");
    return; 
  };
  
  fprintf(fout, "%d\n", out->boundary->num_faces);
  for (i=0; i < out->boundary->num_faces; i++)
  {
    triangle = out->boundary->face + i;
    fprintf(fout, "Tr %d %d %d\n", triangle->a, triangle->b, triangle->c);
  }

  fclose(fout);
  
  // Write all boundary files
  i = 0;
  for (std::set<int>::iterator it = all_markers.begin();
       it != all_markers.end(); it++)
  {
    sprintf (filename, "%s_%s.vtx", basename, boundaries[i]);
    if ((fout=fopen(filename, "wb"))==NULL){
      printf("write error...\n");
      return; 
    };    

    fprintf(fout, "%d\n", out->num_cells);
    fprintf(fout, "intra\n");
    for (std::set<int>::iterator vert_it = vertex_markers[*it].begin();
	 vert_it != vertex_markers[*it].end(); vert_it++)
      fprintf(fout, "%d\n", *vert_it);
    
    // Increase marker index
    i++;
    fclose(fout);
    
  }

}

/*
 * ***************************************************************************
 * Routine:  GemMesh_writeDiffpack
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Write tetrahedralized mesh to disk (DiffPack format)
 * ***************************************************************************
 */
void GemMesh_writeDiffpack(GemMesh* out, char* filename, int num_boundaries, 
			   char* boundaries[])
{
  unsigned int i;
  FILE *fout;
  FETK_VX* vv;
  FETK_SS* ss;
  INT6VECT* hs;

  // Open file
  if ((fout=fopen(filename, "wb"))==NULL){
    printf("write error...\n");
    return; 
  };    

  // Collect face information
  std::vector<std::set<int> > vertex_markers(out->num_vertices);
  std::set<int> all_markers;
  for (i = 0; i < out->num_cells; i++)
  {
    ss = out->ss + i;

    // Count face markers
    if (ss->fa!=0)
    {
      all_markers.insert(ss->fa);
      vertex_markers[ss->nb].insert(ss->fa);
      vertex_markers[ss->nc].insert(ss->fa);
      vertex_markers[ss->nd].insert(ss->fa);
    }
    if (ss->fb!=0)
    {
      all_markers.insert(ss->fb);
      vertex_markers[ss->na].insert(ss->fb);
      vertex_markers[ss->nc].insert(ss->fb);
      vertex_markers[ss->nd].insert(ss->fb);
    }
    if (ss->fc!=0)
    {
      all_markers.insert(ss->fc);
      vertex_markers[ss->na].insert(ss->fc);
      vertex_markers[ss->nb].insert(ss->fc);
      vertex_markers[ss->nd].insert(ss->fc);
    }
    if (ss->fd!=0)
    {
      all_markers.insert(ss->fd);
      vertex_markers[ss->na].insert(ss->fd);
      vertex_markers[ss->nb].insert(ss->fd);
      vertex_markers[ss->nc].insert(ss->fd);
    }
  }

  // Check that number of facet markers equals the number passed boundary names
  if (all_markers.size()!=num_boundaries)
  {
    throw std::runtime_error("Number of passed boundaries does not match number "
			    "of facet markers in mesh."); 
  }

  // Create a map between facet markers in GemMesh and numbered boundary names
  i = 1;
  std::map<int, int> marker_map;
  for (std::set<int>::iterator it = all_markers.begin();
       it != all_markers.end(); it++)
  {
    marker_map[*it] = i++;
  }

  // Output header
  fprintf(fout, "  Finite element mesh (GridFE):\n\n");

  fprintf(fout, "  Number of space dim. =  3  embedded in physical space with dimension 3\n");
  fprintf(fout, "  Number of elements   =  %d\n", out->num_cells);
  fprintf(fout, "  Number of nodes      =  %d\n\n", out->num_vertices);

  fprintf(fout, "  All elements are of the same type : true\n");
  fprintf(fout, "  Max number of nodes in an element: %d\n", out->higher_order ? 10 : 4);

  // FIXME: Make num materials flexible
  fprintf(fout, "  Only one material                : true\n");
  fprintf(fout, "  Lattice data                     ? 0\n\n\n\n");

  // Output boundary indicator names
  fprintf(fout, "  %d boundary indicators:\n", num_boundaries);
  for (i = 0; i < num_boundaries; i++)
  {
    fprintf(fout, "    %s \n", boundaries[i]);
  }
  fprintf(fout, "\n\n");
  
  fprintf(fout, "  Nodal coordinates and nodal boundary indicators,\n"\
                "  the columns contain:\n"\
                "     - node number\n"\
                "     - coordinates\n"\
                "     - no of boundary indicators that are set (ON)\n"\
                "     - the boundary indicators that are set (ON) if any.\n");
  fprintf(fout, "#\n");

  
  // Write all vertices to file
  char buffer [1000];
  for (i = 0; i < out->num_vertices; i++)
  {
    vv = out->vv + i;
    sprintf(buffer, " ");
    for (std::set<int>::iterator it = vertex_markers[i].begin();
	 it != vertex_markers[i].end(); it++)
    {
      sprintf (buffer, "%s %d ", buffer, marker_map[*it]);
    }
    
    // NOTE: Diffpack use 1 as its first vertex number
    fprintf(fout, "%6d (%e, %e, %e)  [%d] %s\n", i+1, vv->x, vv->y, vv->z, \
	    (int)vertex_markers[i].size(), buffer);
  }

  fprintf(fout, "\n\n");
  fprintf(fout, "  Element types and connectivity\n"\
                "  the columns contain:\n"\
                "   - element number\n"\
                "   - element type\n"\
                "   - material number\n"\
                "   - the global node numbers of the nodes in the element.\n");
  fprintf(fout, "#\n");

  // Write all cells to file 
  for (i = 0; i < out->num_cells; i++)
  {
    ss = out->ss + i;

    // Check sign of determinant det(a - b, b - c, c - d)|
    // If negative change or of last two vertices
    //const double ax = out->vv[ss->na].x;
    //const double ay = out->vv[ss->na].y;
    //const double az = out->vv[ss->na].z;
    //
    //const double bx = out->vv[ss->nb].x;
    //const double by = out->vv[ss->nb].y;
    //const double bz = out->vv[ss->nb].z;
    //
    //const double cx = out->vv[ss->nc].x;
    //const double cy = out->vv[ss->nc].y;
    //const double cz = out->vv[ss->nc].z;
    //
    //const double dx = out->vv[ss->nd].x;
    //const double dy = out->vv[ss->nd].y;
    //const double dz = out->vv[ss->nd].z;

    /*
    // Formula for the determinant of tetrahedron from http://mathworld.wolfram.com
    const double det_j = (ax*(by*cz + dy*bz + cy*dz - cy*bz - by*dz - dy*cz) -
		          bx*(ay*cz + dy*az + cy*dz - cy*az - ay*dz - dy*cz) +
		          cx*(ay*bz + dy*az + by*dz - by*az - ay*dz - dy*bz) -
		          dx*(ay*bz + by*cz + cy*az - by*az - cy*bz - ay*cz));

    */
    //const double j11 = bx - ax; 
    //const double j12 = by - ay; 
    //const double j13 = bz - az;
    //const double j21 = cx - ax; 
    //const double j22 = cy - ay; 
    //const double j23 = cz - az;
    //const double j31 = dx - ax; 
    //const double j32 = dy - ay; 
    //const double j33 = dz - az;
    //
    //const double a11 = j22*j33 - j23*j32;  
    //const double a12 = j13*j32 - j12*j33;  
    //const double a13 = j12*j23 - j13*j22;
    //const double a21 = j23*j31 - j21*j33;  
    //const double a22 = j11*j33 - j13*j31;  
    //const double a23 = j13*j21 - j11*j23;
    //const double a31 = j21*j32 - j22*j31;  
    //const double a32 = j12*j31 - j11*j32;  
    //const double a33 = j11*j22 - j12*j21;
    //
    //const double det_j = j11*a11 + j12*a21 + j13*a31;
    //
    //// NOTE: Diffpack use 1 as its first cell number, also note that each 
    //// NOTE: vertex index also need to be bumped with 1
    //if (det_j > 0)
    //  fprintf(fout, "%6d   ElmT4n3D  %3d  %4d %4d %4d %4d \n", i+1, ss->mat, \
    //	      ss->na+1, ss->nb+1, ss->nc+1, ss->nd+1);
    //else
    //  fprintf(fout, "%6d   ElmT4n3D  %3d  %4d %4d %4d %4d \n", i+1, ss->mat, \
    //      ss->na+1, ss->nb+1, ss->nd+1, ss->nc+1);

    if (out->higher_order)
    {
      hs = out->hs + i;
      
      fprintf(fout, "%6d   ElmT10n3D  %3d  %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d \n", 
	      i+1, ss->mat, ss->na+1, ss->nb+1, ss->nd+1, ss->nc+1, hs->a+1, 
	      hs->b+1, hs->c+1, hs->d+1, hs->e+1, hs->f+1);
    }
    else
    {
      fprintf(fout, "%6d   ElmT4n3D  %3d  %4d %4d %4d %4d \n", i+1, ss->mat, \
	      ss->na+1, ss->nb+1, ss->nd+1, ss->nc+1);
    }

  }

  fclose(fout);

}

/*
 * ***************************************************************************
 * Routine:  GemMesh_writeMcsf
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Write tetrahedralized mesh to disk (MCSF formats)
 * ***************************************************************************
 */

void GemMesh_writeMcsf(GemMesh* out, char* filename)
{
  if (out->higher_order)
  {
    printf("MCSF output does not support higher order mesh...\n");
    return; 
  }

  unsigned int i;
  FILE *fout;
  FETK_VX* vv;
  FETK_SS* ss;

  // Open file
  if ((fout=fopen(filename, "wb"))==NULL){
    printf("write error...\n");
    return; 
  };    

  // File head
  fprintf(fout,"mcsf_begin=1; \n \n");
  fprintf(fout,"      dim = %d; \n", out->dim);
  fprintf(fout,"    dimii = %d; \n", out->dimii);
  fprintf(fout," vertices = %d;\n", out->num_vertices);
  fprintf(fout,"simplices = %d;\n", out->num_cells);

  // Node head
  fprintf(fout,"vert=[ \n");
  fprintf(fout,"%%------------------------------------------------------------------------------ \n");
  fprintf(fout,"%%  Node-ID  Chrt        X-Coordinate        Y-coordinate        Z-coordinate    \n");
  fprintf(fout,"%%---------  ----     ----------------     ----------------     ---------------- \n");
  
  // Write all nodes with type to file
  for (i = 0; i < out->num_vertices; i++){
    vv = out->vv + i;
    fprintf(fout, "%10d    %d     %17.10e    %17.10e    %17.10e   \n", i, 
	    vv->chrt, vv->x, vv->y, vv->z);
  }
  
  fprintf(fout,"]; \n");

  // Cell (simplex) head
  fprintf(fout,"simp=[ \n");
  fprintf(fout,"%%------------------------------------------------------------------------------------------ \n");
  fprintf(fout,"%%  Simp-ID Grp    Mat          Face-Types                      Vertex-Numbers \n");
  fprintf(fout,"%%--------- ---    ---    ---------------------  ------------------------------------------- \n");

  for (i = 0; i < out->num_cells; i++) {
    ss = out->ss + i;
    fprintf(fout,"%10d  0  %5d  %5d %5d %5d %5d  %10d %10d %10d %10d \n", i, ss->id, 
	    ss->fa, ss->fb, ss->fc, ss->fd, ss->na, ss->nb, ss->nc, ss->nd);
  }

  fprintf(fout,"]; \n");
  fprintf(fout,"mcsf_end=1; \n");
  fclose(fout);

}

/*
 * ***************************************************************************
 * Routine:  write_mcsf
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Write surface/tetrahedral meshes into disks (MCSF/OFF formats)
 * ***************************************************************************
 */
//void tetgenio::write_mcsf(char* filename, unsigned int* ActiveSite)
//{
//  unsigned int i, j, k;
//  FILE *fout;
//  int *node_attr;
//  long node_index;
//  int *face_type, tetra_node[4];
//  float x,y,z;
//  char file_name[256];
//  int cell_type, interior_node_type, exterior_node_type;
//  float dist;
//  long a,b,c,d;
//  float ax,ay,az;
//  float bx,by,bz;
//  float cx,cy,cz;
//  float dx,dy,dz;
//  int active, exterior_cell_type;
//  long neighbor, neightype;
//  long *map_w_o,*map_w_i,*map_w_t; 
//  
//  printf("Num trifaces: %d\n", numberoftrifaces);
//
//  // Check the type of exterior tetrahedra
//  // In Tetgen, it could be radmonly 1 or 2 (or even other integers if more materials are present)
//  fflush(stdout);
//  exterior_cell_type = -1;
//  for (i = 0; i < numberoftetrahedra; i++) {
//    cell_type = (int)tetrahedronattributelist[i * numberoftetrahedronattributes];
//    for (j = 0; j < numberofcorners; j++) {
//      neighbor = neighborlist[i * 4 + j] - 1;
//      if (neighbor < 0) {
//	exterior_cell_type = cell_type;
//	break;
//      }
//    }
//    if (exterior_cell_type >= 0)
//      break;
//  }
//  printf("exterior_cell_type = %d \n", exterior_cell_type);
//
//  // assign node attributes
//  // Node attributes are assigned as either being situated at the exterior or 
//  // interior domain
//  // For now we go for:
//  // 0 = interior domain
//  // 1 = exterior domain
//  // This might be user specific in future
//  interior_node_type = 0;
//  exterior_node_type = 1;
//
//  node_attr = new int[numberofpoints];
//  for (i = 0; i < numberofpoints; i++) {
//    node_attr[i] = 0;
//  }
//  
//  for (i = 0; i < numberoftetrahedra; i++) {
//    cell_type = (int)tetrahedronattributelist[i * numberoftetrahedronattributes];
//    for (j = 0; j < numberofcorners; j++) {
//      // Note that index in tetrahedronlist starts from 1 !!!!!!!!
//      node_index = tetrahedronlist[i * numberofcorners + j] - 1;
//      
//      // If exterior_cell_type is 1
//      if (exterior_cell_type == 1) {
//	if (cell_type >= 3) cell_type = 2;
//	// FIXME: I do not understand this logic. Slowly getting there... JH
//	node_attr[node_index] = cell_type==1 ? exterior_node_type : interior_node_type;
//      }
//      else if (exterior_cell_type == 2) {
//	if (cell_type >= 3) cell_type = 1;
//	node_attr[node_index] = cell_type==2 ? exterior_node_type : interior_node_type;
//      }
//      else {
//	printf("there are extra materials inside molecular surfaces\n");
//	exit(0);
//      }
//    }
//
//    // check tetra orientation
//    a = tetrahedronlist[i * numberofcorners]-1;
//    b = tetrahedronlist[i * numberofcorners + 1]-1;
//    c = tetrahedronlist[i * numberofcorners + 2]-1;
//    d = tetrahedronlist[i * numberofcorners + 3]-1;
//    ax = pointlist[a * 3];
//    ay = pointlist[a * 3 + 1];
//    az = pointlist[a * 3 + 2];
//    bx = pointlist[b * 3];
//    by = pointlist[b * 3 + 1];
//    bz = pointlist[b * 3 + 2];
//    cx = pointlist[c * 3];
//    cy = pointlist[c * 3 + 1];
//    cz = pointlist[c * 3 + 2];
//    dx = pointlist[d * 3];
//    dy = pointlist[d * 3 + 1];
//    dz = pointlist[d * 3 + 2];
//    bx -= ax;
//    by -= ay;
//    bz -= az;
//    cx -= ax;
//    cy -= ay;
//    cz -= az;
//    dx -= ax;
//    dy -= ay;
//    dz -= az;
//    ax = by*cz-bz*cy;
//    ay = bz*cx-bx*cz;
//    az = bx*cy-by*cx;
//    if (ax*dx+ay*dy+az*dz < 0) {
//      tetrahedronlist[i * numberofcorners + 1] = c;
//      tetrahedronlist[i * numberofcorners + 2] = b;
//      k = neighborlist[i * 4 + 1];
//      neighborlist[i * 4 + 1] = neighborlist[i * 4 + 2];
//      neighborlist[i * 4 + 2] = k;
//    }
//  }
//
//  // output the whole tetra mesh
//  if ((fout=fopen(filename, "wb"))==NULL){
//    printf("write error...\n");
//    exit(0); 
//  };    
//
//  // File head
//  fprintf(fout,"mcsf_begin=1; \n \n");
//  fprintf(fout,"      dim = 3; \n");
//  fprintf(fout,"    dimii = 3; \n");
//  fprintf(fout," vertices = %d;\n", numberofpoints);
//  fprintf(fout,"simplices = %d;\n", numberoftetrahedra);
//
//  // Node head
//  fprintf(fout,"vert=[ \n");
//  fprintf(fout,"%%------------------------------------------------------------------------------ \n");
//  fprintf(fout,"%%  Node-ID  Chrt        X-Coordinate        Y-coordinate        Z-coordinate    \n");
//  fprintf(fout,"%%---------  ----     ----------------     ----------------     ---------------- \n");
//  
//  // Write all nodes with type to file
//  for (i = 0; i < numberofpoints; i++)
//    fprintf(fout, "%10d    %d     %17.10e    %17.10e    %17.10e   \n", i, 
//	    node_attr[i], pointlist[i*3], pointlist[i*3+1], pointlist[i*3+2]);
//  
//  fprintf(fout,"]; \n");
//
//  // Cell (simplex) head
//  fprintf(fout,"simp=[ \n");
//  fprintf(fout,"%%------------------------------------------------------------------------------------------ \n");
//  fprintf(fout,"%%  Simp-ID Grp    Mat          Face-Types                      Vertex-Numbers \n");
//  fprintf(fout,"%%--------- ---    ---    ---------------------  ------------------------------------------- \n");
//
//  face_type  = new int[numberoftetrahedra*4];
//  for (i = 0; i < numberoftetrahedra; i++) {
//    for (j = 0; j < numberofcorners; j++) 
//      // NOTE: Index in tetrahedronlist starts from 1
//      tetra_node[j] = tetrahedronlist[i * numberofcorners + j] - 1;
//    
//    cell_type = (int)tetrahedronattributelist[i * numberoftetrahedronattributes];
//
//    // Check face type
//    for (j = 0; j < numberofcorners; j++) {
//      face_type[4*i+j] = 0; 
//      
//      // NOTE: Index in neighborlist starts from 1
//      neighbor = neighborlist[i * 4 + j] - 1;
//      if (neighbor < 0) {
//	face_type[4*i+j] = 5;
//	if (cell_type != exterior_cell_type)
//	  printf("there might be an error with tet attibutes ....\n");
//      }
//      else {
//	neightype = (int)tetrahedronattributelist[neighbor * numberoftetrahedronattributes];
//	
//	if (ActiveSite != NULL){
//	  if (neightype != cell_type) {
//	    active = 1;
//	    for (k = 0; k < numberofcorners; k++) {
//	      if (k != j && ActiveSite[tetra_node[k]] == 0) 
//		active = 0;
//	    }
//	    if (active)
//	      face_type[4*i+j] = ActiveSite[tetra_node[(j+1)%numberofcorners]];
//	    else 
//	      face_type[4*i+j] = 4;
//	  }
//	}
//      }
//    }
//
//    if (exterior_cell_type == 1) {
//      if (cell_type >= 2) cell_type = 1;
//      else cell_type = 2;
//    }
//
//    else if (exterior_cell_type == 2) {
//      if (cell_type >= 3) cell_type = 2;
//    }
//    else {
//      printf("there are extra materials inside molecular surfaces\n");
//      exit(0);
//    }
//
//    fprintf(fout,"%10d  0  %5d  %5d %5d %5d %5d  %10d %10d %10d %10d \n", i, 0, 
//	    face_type[4*i], face_type[4*i+1], face_type[4*i+2], face_type[4*i+3], 
//	    tetra_node[0], tetra_node[1], tetra_node[2], tetra_node[3]);
//  }
//
//  fprintf(fout,"]; \n");
//  fprintf(fout,"mcsf_end=1; \n");
//  fclose(fout);
//
//  // Release memory
//  delete[] node_attr;
//  delete[] face_type;
//
//}

/*
 * ***************************************************************************
 * Routine:  write_pdb_gem
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Write tetrahedral meshes into "gem_mesh", which can be readily
 *           converted to FETK-supported format
 * ***************************************************************************
 */

GemMesh* GemMesh_fromPdb(tetgenio* out, float radius, float centerx, 
			 float centery, float centerz, char *ActiveSite, 
			 int output_flag)
{
  GemMesh* gem_mesh;
  long i, j, k;
  int *NodeAttr;
  long NodeIndex;
  int *FaceType, TetraNode[4];
  float x,y,z;
  int type;
  float dist;
  long a,b,c,d;
  float ax,ay,az;
  float bx,by,bz;
  float cx,cy,cz;
  float dx,dy,dz;
  int active, ExType;
  long neighbor, neightype;
  long *map_w_o,*map_w_i;
  long out_pnt_num,out_tet_num;
  long in_pnt_num,in_tet_num;


  // initialization and memory allocation
  gem_mesh = (GemMesh *)malloc(sizeof(GemMesh));
  gem_mesh->dim = 3;
  gem_mesh->dimii = 3;
  gem_mesh->vv = (FETK_VX *)malloc(sizeof(FETK_VX)*out->numberofpoints);
  gem_mesh->ss = (FETK_SS *)malloc(sizeof(FETK_SS)*out->numberoftetrahedra);


  // Check the type of exterior tetrahedra
  // In Tetgen, it could be radmonly 1 or 2 (or even other integers if more materials are present)
  ExType = -1;
  for (i = 0; i < out->numberoftetrahedra; i++) {
    type = (int)out->tetrahedronattributelist[i * out->numberoftetrahedronattributes];
    for (j = 0; j < out->numberofcorners; j++) {
      neighbor = out->neighborlist[i * 4 + j] - 1;
      if (neighbor < 0) {
	ExType = type;
	break;
      }
    }
    if (ExType >= 0)
      break;
  }
  printf("ExType = %d \n",ExType);

  // assign node attributes
  NodeAttr = new int[out->numberofpoints];
  for (i = 0; i < out->numberofpoints; i++) {
    if (out->pointmarkerlist[i]) {
      x = out->pointlist[i * 3];
      y = out->pointlist[i * 3 + 1];
      z = out->pointlist[i * 3 + 2];
      dist = sqrt((x-centerx)*(x-centerx)+(y-centery)*(y-centery)+(z-centerz)*(z-centerz));
      if (dist > 0.99*radius)
	NodeAttr[i] = 4;
      else
	NodeAttr[i] = 3;
    }
    else 
      NodeAttr[i] = 0;
  }
  
  for (i = 0; i < out->numberoftetrahedra; i++) {
    type = (int)out->tetrahedronattributelist[i * out->numberoftetrahedronattributes];
    for (j = 0; j < out->numberofcorners; j++) {
      // Note that index in tetrahedronlist starts from 1 !!!!!!!!
      NodeIndex = out->tetrahedronlist[i * out->numberofcorners + j] - 1;
      
      if (NodeAttr[NodeIndex] == 0) {
	if (ExType == 1) {
	  if (type >= 3) type = 2;
	  NodeAttr[NodeIndex] = ((type==1) ? (2):(1));
	}
	else if (ExType == 2) {
	  if (type >= 3) type = 1;
	  NodeAttr[NodeIndex] = type;
	}
	else {
	  printf("there are extra materials inside molecular surfaces\n");
	  return NULL;
	}
      }
    }
    // check tetra orientation
    a = out->tetrahedronlist[i * out->numberofcorners]-1;
    b = out->tetrahedronlist[i * out->numberofcorners + 1]-1;
    c = out->tetrahedronlist[i * out->numberofcorners + 2]-1;
    d = out->tetrahedronlist[i * out->numberofcorners + 3]-1;
    ax = out->pointlist[a * 3];
    ay = out->pointlist[a * 3 + 1];
    az = out->pointlist[a * 3 + 2];
    bx = out->pointlist[b * 3];
    by = out->pointlist[b * 3 + 1];
    bz = out->pointlist[b * 3 + 2];
    cx = out->pointlist[c * 3];
    cy = out->pointlist[c * 3 + 1];
    cz = out->pointlist[c * 3 + 2];
    dx = out->pointlist[d * 3];
    dy = out->pointlist[d * 3 + 1];
    dz = out->pointlist[d * 3 + 2];
    bx -= ax;
    by -= ay;
    bz -= az;
    cx -= ax;
    cy -= ay;
    cz -= az;
    dx -= ax;
    dy -= ay;
    dz -= az;
    ax = by*cz-bz*cy;
    ay = bz*cx-bx*cz;
    az = bx*cy-by*cx;
    if (ax*dx+ay*dy+az*dz < 0) {
      out->tetrahedronlist[i * out->numberofcorners + 1] = c;
      out->tetrahedronlist[i * out->numberofcorners + 2] = b;
      k = out->neighborlist[i * 4 + 1];
      out->neighborlist[i * 4 + 1] = out->neighborlist[i * 4 + 2];
      out->neighborlist[i * 4 + 2] = k;
    }
  }
  
  
  // output the whole tetra mesh  
  if (output_flag == 1) {
    gem_mesh->num_vertices  = out->numberofpoints;
    gem_mesh->num_cells = out->numberoftetrahedra;
  }
  out_pnt_num = 0;
  in_pnt_num = 0;
  for (i = 0; i < out->numberofpoints; i++) {
    if (output_flag == 1) {
      gem_mesh->vv[i].id = i;
      gem_mesh->vv[i].chrt = (NodeAttr[i]==4)?(2):(NodeAttr[i]);
      gem_mesh->vv[i].x = out->pointlist[i * 3];
      gem_mesh->vv[i].y = out->pointlist[i * 3 + 1];
      gem_mesh->vv[i].z = out->pointlist[i * 3 + 2];
    }
    if (NodeAttr[i] != 1) 
      out_pnt_num++;
    if (NodeAttr[i] == 1 || NodeAttr[i] == 3) 
      in_pnt_num++;
  }
  
  out_tet_num = 0;
  in_tet_num = 0;
  FaceType  = new int[out->numberoftetrahedra*4];
  for (i = 0; i < out->numberoftetrahedra; i++) {
    for (j = 0; j < out->numberofcorners; j++) 
      // Note that index in tetrahedronlist starts from 1 !!!!!!!!
      TetraNode[j] = out->tetrahedronlist[i * out->numberofcorners + j] - 1;
    
    type = (int)out->tetrahedronattributelist[i * out->numberoftetrahedronattributes];
    if (type == ExType)
      out_tet_num++;
    else
      in_tet_num++;

    // Check face type
    for (j = 0; j < out->numberofcorners; j++) {
      FaceType[4*i+j] = 0; 
      
      // Note that index in neighborlist starts from 1 !!!!!!!!
      neighbor = out->neighborlist[i * 4 + j] - 1;
      if (neighbor < 0) {
	FaceType[4*i+j] = 5;
	if (type != ExType)
	  printf("there might be an error with tet attibutes ....\n");
      }
      else {
	neightype = (int)out->tetrahedronattributelist[neighbor * out->numberoftetrahedronattributes];
	if (ActiveSite != NULL){
	  if (neightype != type) {
	    active = 1;
	    for (k = 0; k < out->numberofcorners; k++) {
	      if (k != j && ActiveSite[TetraNode[k]] == 0) 
		active = 0;
	    }
	    if (active)
	      FaceType[4*i+j] = ActiveSite[TetraNode[(j+1)%out->numberofcorners]];
	    else 
	      FaceType[4*i+j] = 4;
	  }
	}
      }
    }
    if (ExType == 1) {
      if (type >= 2) type = 1;
      else type = 2;
    }
    else if (ExType == 2) {
      if (type >= 3) type = 2;
    }
    else {
      printf("there are extra materials inside molecular surfaces\n");
      return NULL;
    }
    if (output_flag == 1) {
      gem_mesh->ss[i].id = i;
      gem_mesh->ss[i].grp = 0;
      gem_mesh->ss[i].mat = type;
      gem_mesh->ss[i].fa = FaceType[4*i];
      gem_mesh->ss[i].fb = FaceType[4*i+1];
      gem_mesh->ss[i].fc = FaceType[4*i+2];
      gem_mesh->ss[i].fd = FaceType[4*i+3];
      gem_mesh->ss[i].na = TetraNode[0];
      gem_mesh->ss[i].nb = TetraNode[1];
      gem_mesh->ss[i].nc = TetraNode[2];
      gem_mesh->ss[i].nd = TetraNode[3];
    }
  }


  // output the inner tetra mesh
  if (output_flag == 2) {
    gem_mesh->num_vertices = in_pnt_num;
    gem_mesh->num_cells = in_tet_num;
  }
  j = 0;
  map_w_i = new long[out->numberofpoints];
  for (i = 0; i < out->numberofpoints; i++) {
    map_w_i[i] = j;
    if (NodeAttr[i] == 1 || NodeAttr[i] == 3) {
      if (output_flag == 2) {
	gem_mesh->vv[j].id = j;
	gem_mesh->vv[j].chrt = 0;
	gem_mesh->vv[j].x = out->pointlist[i * 3];
	gem_mesh->vv[j].y = out->pointlist[i * 3 + 1];
	gem_mesh->vv[j].z = out->pointlist[i * 3 + 2];
      }
      j++;
    }
  }
  k = 0;
  for (i = 0; i < out->numberoftetrahedra; i++) {
    type = (int)out->tetrahedronattributelist[i * out->numberoftetrahedronattributes];
    if (type != ExType) {
      for (j = 0; j < out->numberofcorners; j++)
	TetraNode[j] = map_w_i[out->tetrahedronlist[i * out->numberofcorners + j] - 1];
      if (output_flag == 2) {
	gem_mesh->ss[k].id = k;
	gem_mesh->ss[k].grp = 0;
	gem_mesh->ss[k].mat = 0;
	gem_mesh->ss[k].fa = FaceType[i*4];
	gem_mesh->ss[k].fb = FaceType[i*4+1];
	gem_mesh->ss[k].fc = FaceType[i*4+2];
	gem_mesh->ss[k].fd = FaceType[i*4+3];
	gem_mesh->ss[k].na = TetraNode[0];
	gem_mesh->ss[k].nb = TetraNode[1];
	gem_mesh->ss[k].nc = TetraNode[2];
	gem_mesh->ss[k].nd = TetraNode[3]; 
      }
      k++;
    }
  }
  

  // output the outer tetra mesh
  if (output_flag == 3) {
    gem_mesh->num_vertices = out_pnt_num;
    gem_mesh->num_cells = out_tet_num;
  }    
  j = 0;
  map_w_o = new long[out->numberofpoints];
  for (i = 0; i < out->numberofpoints; i++) {
    map_w_o[i] = j;
    if (NodeAttr[i] != 1) {
      if (output_flag == 3) {
	gem_mesh->vv[j].id = j;
	gem_mesh->vv[j].chrt = 0;
	gem_mesh->vv[j].x = out->pointlist[i * 3];
	gem_mesh->vv[j].y = out->pointlist[i * 3 + 1];
	gem_mesh->vv[j].z = out->pointlist[i * 3 + 2];
      }
      j++;
    }
  }
  k = 0;
  for (i = 0; i < out->numberoftetrahedra; i++) {
    type = (int)out->tetrahedronattributelist[i * out->numberoftetrahedronattributes];
    if (type == ExType) {
      for (j = 0; j < out->numberofcorners; j++)
	TetraNode[j] = map_w_o[out->tetrahedronlist[i * out->numberofcorners + j] - 1];
      if (output_flag == 3) {
	gem_mesh->ss[k].id = k;
	gem_mesh->ss[k].grp = 0;
	gem_mesh->ss[k].mat = 0;
	gem_mesh->ss[k].fa = FaceType[i*4];
	gem_mesh->ss[k].fb = FaceType[i*4+1];
	gem_mesh->ss[k].fc = FaceType[i*4+2];
	gem_mesh->ss[k].fd = FaceType[i*4+3];
	gem_mesh->ss[k].na = TetraNode[0];
	gem_mesh->ss[k].nb = TetraNode[1];
	gem_mesh->ss[k].nc = TetraNode[2];
	gem_mesh->ss[k].nd = TetraNode[3]; 
      }
      k++;
    }
  }
  
  free(NodeAttr);
  free(FaceType);
  free(map_w_i);
  free(map_w_o);

  return gem_mesh;
}


/*
 * ***************************************************************************
 * Routine:  write_pdb_mcsf
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Write surface/tetrahedral meshes into disks (MCSF/OFF formats)
 * ***************************************************************************
 */
// FIXME: Depricated. Use write to GemMesh and then to file
//void tetgenio::write_pdb_mcsf(char *inputfilename, float radius, float centerx, 
//			  float centery, float centerz, char *ActiveSite, 
//			  int output_flag)
//{
//  long i, j, k;
//  FILE *fout;
//  int *NodeAttr;
//  long NodeIndex;
//  int *FaceType, TetraNode[4];
//  float x,y,z;
//  char file_name[256];
//  int type;
//  float dist;
//  long a,b,c,d;
//  float ax,ay,az;
//  float bx,by,bz;
//  float cx,cy,cz;
//  float dx,dy,dz;
//  int active, ExType;
//  long neighbor, neightype;
//  long out_pnt_num,out_tet_num;
//  long in_pnt_num,in_tet_num;
//  long bem_pnt_num,bem_tri_num;
//  long *map_w_o,*map_w_i,*map_w_t; 
//  
//
//  
//  // Check the type of exterior tetrahedra
//  // In Tetgen, it could be radmonly 1 or 2 (or even other integers if more materials are present)
//  ExType = -1;
//  for (i = 0; i < numberoftetrahedra; i++) {
//    type = (int)tetrahedronattributelist[i * numberoftetrahedronattributes];
//    for (j = 0; j < numberofcorners; j++) {
//      neighbor = neighborlist[i * 4 + j] - 1;
//      if (neighbor < 0) {
//	ExType = type;
//	break;
//      }
//    }
//    if (ExType >= 0)
//      break;
//  }
//  printf("ExType = %d \n",ExType);
//
//  // assign node attributes
//  bem_pnt_num = 0;
//  NodeAttr = new int[numberofpoints];
//  for (i = 0; i < numberofpoints; i++) {
//    if (pointmarkerlist[i]) {
//      x = pointlist[i * 3];
//      y = pointlist[i * 3 + 1];
//      z = pointlist[i * 3 + 2];
//      dist = sqrt((x-centerx)*(x-centerx)+(y-centery)*(y-centery)+(z-centerz)*(z-centerz));
//      if (dist > 0.99*radius)
//	NodeAttr[i] = 4;
//      else {
//	NodeAttr[i] = 3;
//	bem_pnt_num++;
//      }
//    }
//    else 
//      NodeAttr[i] = 0;
//  }
//  
//  for (i = 0; i < numberoftetrahedra; i++) {
//    type = (int)tetrahedronattributelist[i * numberoftetrahedronattributes];
//    for (j = 0; j < numberofcorners; j++) {
//      // Note that index in tetrahedronlist starts from 1 !!!!!!!!
//      NodeIndex = tetrahedronlist[i * numberofcorners + j] - 1;
//      
//      if (NodeAttr[NodeIndex] == 0) {
//	if (ExType == 1) {
//	  if (type >= 3) type = 2;
//	  NodeAttr[NodeIndex] = ((type==1) ? (2):(1));
//	}
//	else if (ExType == 2) {
//	  if (type >= 3) type = 1;
//	  NodeAttr[NodeIndex] = type;
//	}
//	else {
//	  printf("there are extra materials inside molecular surfaces\n");
//	  exit(0);
//	}
//      }
//    }
//    // check tetra orientation
//    a = tetrahedronlist[i * numberofcorners]-1;
//    b = tetrahedronlist[i * numberofcorners + 1]-1;
//    c = tetrahedronlist[i * numberofcorners + 2]-1;
//    d = tetrahedronlist[i * numberofcorners + 3]-1;
//    ax = pointlist[a * 3];
//    ay = pointlist[a * 3 + 1];
//    az = pointlist[a * 3 + 2];
//    bx = pointlist[b * 3];
//    by = pointlist[b * 3 + 1];
//    bz = pointlist[b * 3 + 2];
//    cx = pointlist[c * 3];
//    cy = pointlist[c * 3 + 1];
//    cz = pointlist[c * 3 + 2];
//    dx = pointlist[d * 3];
//    dy = pointlist[d * 3 + 1];
//    dz = pointlist[d * 3 + 2];
//    bx -= ax;
//    by -= ay;
//    bz -= az;
//    cx -= ax;
//    cy -= ay;
//    cz -= az;
//    dx -= ax;
//    dy -= ay;
//    dz -= az;
//    ax = by*cz-bz*cy;
//    ay = bz*cx-bx*cz;
//    az = bx*cy-by*cx;
//    if (ax*dx+ay*dy+az*dz < 0) {
//      tetrahedronlist[i * numberofcorners + 1] = c;
//      tetrahedronlist[i * numberofcorners + 2] = b;
//      k = neighborlist[i * 4 + 1];
//      neighborlist[i * 4 + 1] = neighborlist[i * 4 + 2];
//      neighborlist[i * 4 + 2] = k;
//    }
//  }
//  
//
//  // output the whole tetra mesh  
//  if (output_flag == 1) {
//    sprintf(file_name, "%s.output.all.m", inputfilename);
//    if ((fout=fopen(file_name, "wb"))==NULL){
//      printf("write error...\n");
//      exit(0); 
//    };    
//    fprintf(fout,"mcsf_begin=1; \n \n");
//    fprintf(fout,"      dim = 3; \n");
//    fprintf(fout,"    dimii = 3; \n");
//    fprintf(fout," vertices = %d;\n", numberofpoints);
//    fprintf(fout,"simplices = %d;\n", numberoftetrahedra);
//    fprintf(fout,"vert=[ \n");
//    fprintf(fout,"%%------------------------------------------------------------------------------ \n");
//    fprintf(fout,"%%  Node-ID  Chrt        X-Coordinate        Y-coordinate        Z-coordinate    \n");
//    fprintf(fout,"%%---------  ----     ----------------     ----------------     ---------------- \n");
//  }
//  
//  out_pnt_num = 0;
//  in_pnt_num = 0;
//  for (i = 0; i < numberofpoints; i++) {
//    if (output_flag == 1)
//      fprintf(fout, "%10ld    %d     %17.10e    %17.10e    %17.10e   \n",i,((NodeAttr[i]==4)?(2):(NodeAttr[i])),
//	      pointlist[i * 3],pointlist[i * 3 + 1],pointlist[i * 3 + 2]);
//    if (NodeAttr[i] != 1) 
//      out_pnt_num++;
//    if (NodeAttr[i] == 1 || NodeAttr[i] == 3) 
//      in_pnt_num++;
//  }
//  
//  if (output_flag == 1) {
//    fprintf(fout,"]; \n");
//    fprintf(fout,"simp=[ \n");
//    fprintf(fout,"%%------------------------------------------------------------------------------------------ \n");
//    fprintf(fout,"%%  Simp-ID Grp    Mat          Face-Types                      Vertex-Numbers \n");
//    fprintf(fout,"%%--------- ---    ---    ---------------------  ------------------------------------------- \n");
//  }
//
//  out_tet_num = 0;
//  in_tet_num = 0;
//  bem_tri_num = 0;
//  FaceType  = new int[numberoftetrahedra*4];
//  for (i = 0; i < numberoftetrahedra; i++) {
//    for (j = 0; j < numberofcorners; j++) 
//      // Note that index in tetrahedronlist starts from 1 !!!!!!!!
//      TetraNode[j] = tetrahedronlist[i * numberofcorners + j] - 1;
//    
//    type = (int)tetrahedronattributelist[i * numberoftetrahedronattributes];
//    if (type == ExType)
//      out_tet_num++;
//    else
//      in_tet_num++;
//
//    // Check face type
//    for (j = 0; j < numberofcorners; j++) {
//      FaceType[4*i+j] = 0; 
//      
//      // Note that index in neighborlist starts from 1 !!!!!!!!
//      neighbor = neighborlist[i * 4 + j] - 1;
//      if (neighbor < 0) {
//	FaceType[4*i+j] = 5;
//	if (type != ExType)
//	  printf("there might be an error with tet attibutes ....\n");
//      }
//      else {
//	neightype = (int)tetrahedronattributelist[neighbor * numberoftetrahedronattributes];
//	
//	if (ActiveSite != NULL){
//	  if (neightype != type) {
//	    active = 1;
//	    for (k = 0; k < numberofcorners; k++) {
//	      if (k != j && ActiveSite[TetraNode[k]] == 0) 
//		active = 0;
//	    }
//	    if (active)
//	      FaceType[4*i+j] = ActiveSite[TetraNode[(j+1)%numberofcorners]];
//	    else 
//	      FaceType[4*i+j] = 4;
//	  }
//	}
//      }
//      if (type == ExType && (FaceType[4*i+j] != 0 && FaceType[4*i+j] != 5))
//	bem_tri_num++;
//      
//    }
//    if (ExType == 1) {
//      if (type >= 2) type = 1;
//      else type = 2;
//    }
//    else if (ExType == 2) {
//      if (type >= 3) type = 2;
//    }
//    else {
//      printf("there are extra materials inside molecular surfaces\n");
//      exit(0);
//    }
//    if (output_flag == 1)
//      fprintf(fout,"%10ld  0  %5d  %5d %5d %5d %5d  %10d %10d %10d %10d \n",i,type,((FaceType[4*i]==5)?(5):(0)),
//	      ((FaceType[4*i+1]==5)?(5):(0)),((FaceType[4*i+2]==5)?(5):(0)),((FaceType[4*i+3]==5)?(5):(0)),
//	      TetraNode[0],TetraNode[1],TetraNode[2],TetraNode[3]);
//  }
//  if (output_flag == 1) {
//    fprintf(fout,"]; \n");
//    fprintf(fout,"mcsf_end=1; \n");
//    fclose(fout);
//  }
//
//
//  // output the inner tetra mesh
//  if (output_flag == 2) {
//    sprintf(file_name, "%s.output.in.m", inputfilename);
//    if ((fout=fopen(file_name, "wb"))==NULL){
//      printf("write error...\n");
//      exit(0);
//    };
//    fprintf(fout,"mcsf_begin=1; \n \n");
//    fprintf(fout,"      dim = 3; \n");
//    fprintf(fout,"    dimii = 3; \n");
//    fprintf(fout," vertices = %ld;\n", in_pnt_num);
//    fprintf(fout,"simplices = %ld;\n", in_tet_num);
//    fprintf(fout,"vert=[ \n");
//    fprintf(fout,"%%------------------------------------------------------------------------------------ \n");
//    fprintf(fout,"%%  Node-ID        Chrt         X-Coordinate         Y-coordinate         Z-coordinate   \n");
//    fprintf(fout,"%%---------  ----------       ----------------   ----------------     ---------------- \n");
//  }
//  j = 0;
//  map_w_i = new long[numberofpoints];
//  for (i = 0; i < numberofpoints; i++) {
//    map_w_i[i] = j;
//    if (NodeAttr[i] == 1 || NodeAttr[i] == 3) {
//      if (output_flag == 2)
//	fprintf(fout, "%10ld  %10d     %17.10e    %17.10e    %17.10e   \n",j,0,
//		pointlist[i * 3],pointlist[i * 3 + 1],pointlist[i * 3 + 2]);
//      j++;
//    }
//  }
//  if (output_flag == 2) {
//    fprintf(fout,"]; \n");
//    fprintf(fout,"simp=[ \n");
//    fprintf(fout,"%%--------------------------------------------------------------------------------------- \n");
//    fprintf(fout,"%%  Simp-ID Grp Mat          Face-Types                      Vertex-Numbers \n");
//    fprintf(fout,"%%--------- --- ---    ---------------------  ------------------------------------------- \n");
//  }
//  k = 0;
//  for (i = 0; i < numberoftetrahedra; i++) {
//    type = (int)tetrahedronattributelist[i * numberoftetrahedronattributes];
//    if (type != ExType) {
//      for (j = 0; j < numberofcorners; j++)
//	TetraNode[j] = map_w_i[tetrahedronlist[i * numberofcorners + j] - 1];
//      if (output_flag == 2)
//	fprintf(fout,"%10ld  0   0  %5d %5d %5d %5d  %10d %10d %10d %10d \n",k,FaceType[i*4],FaceType[i*4+1],
//		FaceType[i*4+2],FaceType[i*4+3],TetraNode[0],TetraNode[1],TetraNode[2],TetraNode[3]);
//      k++;
//    }
//  }
//  if (output_flag == 2) {
//    fprintf(fout,"]; \n");
//    fprintf(fout,"mcsf_end=1; \n");
//    fclose(fout);
//  }
//  
//
//  // output the outer tetra mesh
//  if (output_flag == 3) {
//    sprintf(file_name, "%s.output.out.m", inputfilename);
//    if ((fout=fopen(file_name, "wb"))==NULL){
//      printf("write error...\n");
//      exit(0);
//    };
//    fprintf(fout,"mcsf_begin=1; \n \n");
//    fprintf(fout,"      dim = 3; \n");
//    fprintf(fout,"    dimii = 3; \n");
//    fprintf(fout," vertices = %ld;\n", out_pnt_num);
//    fprintf(fout,"simplices = %ld;\n", out_tet_num);
//    fprintf(fout,"vert=[ \n");
//    fprintf(fout,"%%------------------------------------------------------------------------------------ \n");
//    fprintf(fout,"%%  Node-ID        Chrt         X-Coordinate         Y-coordinate         Z-coordinate   \n");
//    fprintf(fout,"%%---------  ----------       ----------------   ----------------     ---------------- \n");
//  }    
//  j = 0;
//  map_w_o = new long[numberofpoints];
//  for (i = 0; i < numberofpoints; i++) {
//    map_w_o[i] = j;
//    if (NodeAttr[i] != 1) {
//      if (output_flag == 3)
//	fprintf(fout, "%10ld  %10d     %17.10e    %17.10e    %17.10e   \n",j,0,
//		pointlist[i * 3],pointlist[i * 3 + 1],pointlist[i * 3 + 2]);
//      j++;
//    }
//  }
//  if (output_flag == 3) {
//    fprintf(fout,"]; \n");
//    fprintf(fout,"simp=[ \n");
//    fprintf(fout,"%%--------------------------------------------------------------------------------------- \n");
//    fprintf(fout,"%%  Simp-ID Grp Mat          Face-Types                      Vertex-Numbers \n");
//    fprintf(fout,"%%--------- --- ---    ---------------------  ------------------------------------------- \n");
//  }
//  k = 0;
//  for (i = 0; i < numberoftetrahedra; i++) {
//    type = (int)tetrahedronattributelist[i * numberoftetrahedronattributes];
//    if (type == ExType) {
//      for (j = 0; j < numberofcorners; j++)
//	TetraNode[j] = map_w_o[tetrahedronlist[i * numberofcorners + j] - 1];
//      if (output_flag == 3)
//	fprintf(fout,"%10ld  0   0  %5d %5d %5d %5d  %10d %10d %10d %10d \n",k,FaceType[i*4],FaceType[i*4+1],
//		FaceType[i*4+2],FaceType[i*4+3],TetraNode[0],TetraNode[1],TetraNode[2],TetraNode[3]);
//      k++;
//    }
//  }
//  if (output_flag == 3) {
//    fprintf(fout,"]; \n");
//    fprintf(fout,"mcsf_end=1; \n");
//    fclose(fout);
//  }



  // the following may be useful for boundary element analysis 
  /*
  // output the mapping file from outer domain to whole domain
  sprintf(file_name, "%s.output.map.m", inputfilename);
  if ((fout=fopen(file_name, "wb"))==NULL){
    printf("write error...\n");
    exit(0);
  };
  j = 0;
  for (i = 0; i < numberofpoints; i++) {
    if (NodeAttr[i] != 1) {
      fprintf(fout,"%10d %10d \n",j+1,i); 
      j++;
    }
  }
  fclose(fout);
  
  
  // output the molecular surface mesh in .m
  sprintf(file_name, "%s.output.bem.m", inputfilename);
  if ((fout=fopen(file_name, "wb"))==NULL){
    printf("write error...\n");
    exit(0);
  };
  fprintf(fout,"%10d %10d 0 \n",bem_pnt_num,bem_tri_num);
  map_w_t = new long[numberofpoints];
  j = 0;
  for (i = 0; i < numberofpoints; i++) {
    map_w_t[i] = j;
    if (NodeAttr[i] == 3) {
      fprintf(fout,"%13.4f %13.4f %13.4f %10d \n",pointlist[i * 3],pointlist[i * 3 + 1],
	      pointlist[i * 3 + 2],map_w_o[i]+1);
      j++;
    }
  }
  for (i = 0; i < numberoftetrahedra; i++) {
    type = (int)tetrahedronattributelist[i * numberoftetrahedronattributes];
    if (type == ExType) {
      for (j = 0; j < numberofcorners; j++) {
	if (FaceType[4*i+j] != 0 && FaceType[4*i+j] != 5) {
	  a = 0;
	  for (k = 0; k < numberofcorners; k++) {
	    if (k != j) { 
	      TetraNode[a] = map_w_t[tetrahedronlist[i * numberofcorners + k] - 1];
	      a++;
	    }
	  }
	  if (j == 0 || j == 2)
	    fprintf(fout,"%10d %10d %10d \n",TetraNode[0]+1,TetraNode[1]+1,TetraNode[2]+1);
	  else if (j == 1 || j == 3)
	    fprintf(fout,"%10d %10d %10d \n",TetraNode[0]+1,TetraNode[2]+1,TetraNode[1]+1);
	}
      }
    }
  }
  fclose(fout);
  
  
  // randomly assign colors to different types of active sites
  float red[256],green[256],blue[256];
  for (i = 0; i < 256; i++) {
    red[i] = (float) rand()/(float)RAND_MAX;
    green[i] = (float) rand()/(float)RAND_MAX;
    blue[i] = (float) rand()/(float)RAND_MAX;
  } 
  // non-active site color
  red[4] = 0.5;
  green[4] = 1;
  blue[4] = 1;
  // active site colors
  red[1] = 1;
  green[1] = 1;
  blue[1] = 0.1;
  red[3] = 0.1;
  green[3] = 0.1;
  blue[3] =  1;
  red[7] = 0.1;
  green[7] = 1;
  blue[7] = 0.1;
  red[9] = 1;
  green[9] = 0.1;
  blue[9] = 0.1;
  red[11] = 1;
  green[11] = 0.1;
  blue[11] = 1;
  red[13] = 1;
  green[13] = 0.5;
  blue[13] = 1;
  // output the molecular surface mesh in .off
  sprintf(file_name, "%s.output.bem.off", inputfilename);
  if ((fout=fopen(file_name, "wb"))==NULL){
    printf("write error...\n");
    exit(0);
  };
  fprintf(fout,"OFF\n%10d %10d 0 \n",bem_pnt_num,bem_tri_num);
  map_w_t = new long[numberofpoints];
  j = 0;
  for (i = 0; i < numberofpoints; i++) {
    map_w_t[i] = j;
    if (NodeAttr[i] == 3) {
      fprintf(fout,"%13.4f %13.4f %13.4f \n",pointlist[i * 3],pointlist[i * 3 + 1],
	      pointlist[i * 3 + 2]);
      j++;
    }
  }
  for (i = 0; i < numberoftetrahedra; i++) {
    type = (int)tetrahedronattributelist[i * numberoftetrahedronattributes];
    if (type == ExType) {
      for (j = 0; j < numberofcorners; j++) {
	if (FaceType[4*i+j] != 0 && FaceType[4*i+j] != 5) { 
	  a = 0;
	  for (k = 0; k < numberofcorners; k++) {
	    if (k != j) { 
	      TetraNode[a] = map_w_t[tetrahedronlist[i * numberofcorners + k] - 1];
	      a++;
	    }
	  }
	  if (j == 0 || j == 2)
	    fprintf(fout,"3 %10d %10d %10d %f %f %f \n",TetraNode[0],TetraNode[1],TetraNode[2],
		    red[FaceType[4*i+j]],green[FaceType[4*i+j]],blue[FaceType[4*i+j]]);
	  else if (j == 1 || j == 3)
	    fprintf(fout,"3 %10d %10d %10d %f %f %f \n",TetraNode[0],TetraNode[2],TetraNode[1],
		    red[FaceType[4*i+j]],green[FaceType[4*i+j]],blue[FaceType[4*i+j]]);
	}
      }
    }
  }
  fclose(fout);
  */
//}

/*
 * ***************************************************************************
 * Routine:  ReadActiveSiteFile
 *
 * Author:   Johan Hake (joha.hake@gmail.com), Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Read a file containing information about active sites
 * ***************************************************************************
 */
void ReadActiveSiteFile(char* input_site, 
			unsigned int& num_spheres, 
			ATOM*& sphere_list, unsigned int*& sphere_markers)
//  ,float*& area_constraint_list)
{
  FILE *fout;
  char buf[1024];
  int i;
  unsigned int marker, c;
  float x, y, z, radius, area_constraint;
  char* ret_gets;
  int ret_scan;
  
  printf("Reading active site file\n");
  if ((fout=fopen(input_site, "r"))==NULL){
    printf("read error...\n");
    return;
  };
  while ((c=fgetc(fout)) == '#')
    ret_gets = fgets(buf, 1024, fout);
  ungetc(c, fout);
  ret_scan = fscanf(fout,"%d\n", &num_spheres);
  sphere_list = (ATOM*)malloc(sizeof(ATOM)*num_spheres);
  sphere_markers  = (unsigned int*)malloc(sizeof(unsigned int)*num_spheres);
  //area_constraint_list = (float*)malloc(sizeof(float)*num_spheres);
  for (i = 0; i < num_spheres; i++) {
    //fscanf(fout, "%f %f %f %f %d\n", &x, &y, &z, &radius, &marker, &area_constraint);
    ret_scan = fscanf(fout, "%f %f %f %f %d\n", &x, &y, &z, &radius, &marker);
    sphere_list[i].x = x;
    sphere_list[i].y = y;
    sphere_list[i].z = z;
    sphere_list[i].radius = radius;
    sphere_markers[i] = marker;
    //area_constraint_list[i] = area_constraint;
  }
  fclose(fout);
}

