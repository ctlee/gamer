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
 * File:     GenerateMesh.C    < ... >
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  The 'main' function of GenerateMesh
 * ***************************************************************************
 */

#include <gamer/gamer.h>

const char* supported_formats[] = {"mcsf", "dolfin"};
const char default_mesh_format[] = "mcsf";
const unsigned int num_supported_formats = 2;
const char tetgen_default_params[] = "nnepAAYM";
const char tetgen_default_quality_params[] = "qq20";

/*
 * ***************************************************************************
 * Routine:      < ... >
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Assign active sites from a list of user-defined spheres
 * ***************************************************************************
 */
//void SetAreaConstraint(SurfaceMesh* surfmesh, ATOM* sphere_list, 
//		       unsigned int* marker_list, int num_spheres,
//		       float* area_constraint_list,
//		       int& num_facet_area_constraints,
//		       double* facet_area_constraints)
//{
//  unsigned int v, sphere, face, a, b, c, i;
//  unsigned int* vertex_spheres, *vertex_markers;
//
//  float x, y, z;
//  float center_x, center_y, center_z;
//  float dist, radius;
//
//  double* tmp_face_area_constraints;
//  unsigned char face_marker;
//
//  // Initalize vertex marker and face marker list memory
//  vertex_spheres = (unsigned int*)malloc(sizeof(unsigned int)*surfmesh->nv);
//  vertex_markers = (unsigned int*)malloc(sizeof(unsigned int)*surfmesh->nv);
//  tmp_face_area_constraints = (double*)malloc(sizeof(double)*surfmesh->nf*2);
//
//  // Reset data
//  for (v = 0; v < surfmesh->nv; v++)
//  {
//    vertex_spheres[v] = 0;
//    vertex_markers[v] = 0;
//  }
//
//  num_facet_area_constraints = 0;
//
//  // Walk through all spheres and mark the vertices
//  for (sphere = 0; sphere < num_spheres; sphere++) {
//    center_x = sphere_list[sphere].x;
//    center_y = sphere_list[sphere].y;
//    center_z = sphere_list[sphere].z;
//    radius   = sphere_list[sphere].radius;
//
//    for (v = 0; v < surfmesh->nv; v++) {
//      x = surfmesh->vertex[v].x;
//      y = surfmesh->vertex[v].y;
//      z = surfmesh->vertex[v].z;
//      dist = sqrt((x-center_x)*(x-center_x)+
//		  (y-center_y)*(y-center_y)+
//		  (z-center_z)*(z-center_z));
//      if (dist < radius){
//	vertex_markers[v] = marker_list[sphere];
//	vertex_spheres[v] = sphere;
//      }
//    }
//  }
//
//  // Walk through the facets and collect area constraint info
//  for (face = 0; face < surfmesh->nf; face++){
//
//    // Get the vertices
//    a = surfmesh->face[face].a;
//    b = surfmesh->face[face].b;
//    c = surfmesh->face[face].c;
//
//    // Get the marker of the vertices and check that all are the same
//    face_marker = vertex_markers[a];
//    if (face_marker != 0 || face_marker != vertex_markers[b])
//      continue;
//    if (face_marker != vertex_markers[c])
//      continue;
//    
//    // If all vertices on a face have the same marker, 
//    // then is the area contraint applied to that face
//    tmp_face_area_constraints[num_facet_area_constraints*2] = face;
//    tmp_face_area_constraints[num_facet_area_constraints*2+1] = \
//      area_constraint_list[vertex_spheres[a]];
//
//    num_facet_area_constraints++;
//  }
//
//  // Create return array
//  facet_area_constraints = new double[num_facet_area_constraints*2];
//
//  // Copy tmp data to return array
//  for (i = 0; i < num_facet_area_constraints*2; i++)
//    facet_area_constraints[i] = tmp_face_area_constraints[i];
//
//  // Free memory
//  free(vertex_markers);
//  free(tmp_face_area_constraints);
//
//}

/*
 * ***************************************************************************
 * Routine:  GenerateMeshFromSurfaceMesh    < ... >
 *
 * Author:   Johan Hake (hake.dev@gmail.com), Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Tetrahedralize a surface mesh
 * ***************************************************************************
 */
void GenerateMeshFromSurfaceMesh(char* input_name, char* input_site, 
				 char* mesh_format, char* tetgen_params)
{
  FILE *fout;

  time_t t0, t1, t2;
  SurfaceMesh* surfmesh = NULL;
  char filename[256], basename[256];
  char* basename_end_pointer;
  char* suffix_end_pointer;
  unsigned int basename_index;

  unsigned int num_spheres = 0;
  ATOM *sphere_list = NULL;
  unsigned int* marker_list = NULL;
  unsigned char* active_sites = NULL;
  float* area_constraint_list = NULL;
  int dummy_count = -1;
  bool use_facemarkers = false;

  tetgenio in, out;
  tetgenio::facet *f;
  tetgenio::polygon *p;

  unsigned int i, j;

  // Read any active sites from file
  if (input_site != NULL)
    ReadActiveSiteFile(input_site, num_spheres, sphere_list, marker_list);
		       //   area_constraint_list);
  
  // Read surface mesh file together with basename
  strcpy(filename, input_name);

  // Check if input is a .off file
  basename_end_pointer = strstr(filename, ".off");
  if (basename_end_pointer == NULL){
    printf("Provide a mesh in \'.off\' format.\n");
    exit(1);
  }
  
  // Get the basename
  basename_index = basename_end_pointer-filename;
  strncpy (basename, filename, basename_index);
  basename[basename_index] = '\0';

  // Load surface meshes in OFF format
  surfmesh = SurfaceMesh_readOFF(input_name);

  printf("Begin generating tetrahedral mesh.\n");
  (void)time(&t0);

  // Convert to a generalized tetrahedral mesh structure
  GemMesh* gem_out = GemMesh_fromSurfaceMesh(surfmesh, tetgen_params);//GemMesh_fromTetgen(out);

  // Save output mesh in given format
  if (strcmp(mesh_format, "mcsf")==0)
  {
    fflush(stdout);
    sprintf(filename, "%s.m", basename);
    printf("Writing results to: %s\n", filename);
    GemMesh_writeMcsf(gem_out, filename);
  }
  
  // Save output mesh in given format
  else if (strcmp(mesh_format, "dolfin")==0)
  {
    fflush(stdout);
    sprintf(filename, "%s.xml", basename);
    printf("Writing results to: %s\n", filename);
    GemMesh_writeDolfin(gem_out, filename);
  }

  // Release memory
  GemMesh_dtor(gem_out);
  SurfaceMesh_dtor(surfmesh);
}


bool check_format(char* format)
{
  for (unsigned int i=0; i<num_supported_formats; i++)
    if (strcmp(format, supported_formats[i]) == 0)
      return true;
  return false;
}

void get_supported_formats(char* formats)
{
  for (unsigned int i=0; i<num_supported_formats; i++)
    if (i == 0)
      sprintf(formats, "\'%s\'", supported_formats[i]);
    else if (i < num_supported_formats-1)
      sprintf(formats, "%s, \'%s\'", formats, supported_formats[i]);
    else
      sprintf(formats, "%s and \'%s\'", formats, supported_formats[i]);
}

void long_usage()
{
  char mesh_formats[1024];
  get_supported_formats(mesh_formats);
  printf("Usage: GenerateMesh [OPTION] FILE\n");
  printf("Generates a tetrahedral mesh from a surface mesh given in FILE");
  printf("\n");
  printf(" FILE - A surface mesh given in .off format\n");
  printf("\n");
  printf("Options:\n");
  printf(" --tetgen-params PAR  Quality parameters which are passed to tetgen:\n");
  printf("                      Default: \'%s\'\n", tetgen_default_quality_params);
  printf("                      Example: \'q2a0.5\'\n");
  printf("                      In addition are these parameters automatically used:\n");
  printf("                      \'%s\'\n", tetgen_default_params);
  printf(" --format FORMAT      The output mesh format. Supported formats:\n");
  printf("                      %s\n", mesh_formats);
  printf("                      Default: \'%s\'\n", default_mesh_format);
  printf(" --active-sites FILE  A file with active sites which will not be coarsed\n" 
	 "                      during quality improvement. \n"
	 "                      See README for specifications.\n"); 
  printf(" --help               Print this help.\n");
}

void option_error()
{
  printf("Error in parsing command line:\nSee \'GenerateMesh --help\' for more information\n");
  exit(1);
}

void parse_comand_line(int argc, char* argv[], char*& input_filename, 
		       char*& input_sites, char*& mesh_format, char*& tetgen_params)
{
  int i = 1;

  if (argc == 1)
  {
    printf("Usage: GenerateMesh [OPTION] FILE\n");
    exit(0);
  }

  // Parse command line
  while (i < argc)
  {
    // Parse --help
    if (strcmp (argv[i], "--help") == 0)
    {
      long_usage();
      exit(0);
    }

    // Parse --format option
    if (strcmp (argv[i], "--format") == 0)
    {
      i++;
      if (i == argc)
	option_error();
      mesh_format = argv[i];
      if (!check_format(mesh_format))
	option_error();
      i++;
      if (i == argc)
	option_error();
      
      continue;
    }

    // Parse --format option
    if (strcmp (argv[i], "--tetgen-params") == 0)
    {
      i++;
      if (i == argc)
	option_error();
      sprintf(tetgen_params, "%s%s", tetgen_default_params, argv[i]);
      i++;
      if (i == argc)
	option_error();
      
      continue;
    }

    // Parse --active-sites option
    if (strcmp (argv[i], "--active-sites") == 0)
    {
      i++;
      if (i == argc)
	option_error();
      input_sites = argv[i];
      i++;
      if (i == argc)
	option_error();
      
      continue;
    }
    
    // Parse filename
    input_filename = argv[i];
    i++;
    if (i != argc)
      option_error();
    
  }

}

int main(int argc, char *argv[])
{
  char* input_sites = NULL;
  char* input_filename = NULL;
  char* mesh_format = NULL;
  char* passed_mesh_format = (char*)malloc(sizeof(char)*16);;
  char* tetgen_params = (char*)malloc(sizeof(char)*16);

  // Parse command line
  parse_comand_line(argc, argv, input_filename, input_sites, 
		    mesh_format, tetgen_params);

  // Set default tetgen params
  sprintf(tetgen_params, "%s%s", tetgen_default_params, tetgen_default_quality_params);
  if (mesh_format == NULL)
    sprintf(passed_mesh_format, "%s", default_mesh_format);
  else
    sprintf(passed_mesh_format, "%s", mesh_format);

  // Call GenerateMesh routine
  GenerateMeshFromSurfaceMesh(input_filename, input_sites, 
			      passed_mesh_format, tetgen_params);
  
  // Free memory
  free(tetgen_params);
  free(passed_mesh_format);

  return(0);
}
