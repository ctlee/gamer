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
 * File:     ImproveSurfMesh.C
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  The 'main' function of ImproveSurfMesh
 * ***************************************************************************
 */


#include <gamer/gamer.h>

// Constants and default values
const char OUT_SUFFIX[] = "_improved_";
const unsigned char OUT_SUFFIX_LENGTH = 10;
const unsigned int MAX_ITER = 6;
const unsigned int NUM_ITER = 1;
const unsigned int MAX_MIN_ANGLE = 15;
const unsigned int MIN_MAX_ANGLE = 150;
const unsigned int PRESERVE_RIDGES = 0;
const float FLAT_RATE = 0.016;
const float DENSE_RATE = 1.5;
const float MAX_MAX_NORMAL_ANGLE = -1.0;

// flag for changed mesh
bool mesh_changed = false;

// Action enums
enum ImproveAction {no_action, refine, smooth, normal_smooth, coarse_dense, 
		    coarse_flat, correct_normals, stats};

// Action struct
struct ActionItem {
  ImproveAction action;
  unsigned int max_iter;
  unsigned int num_iter;
  unsigned int max_min_angle;
  unsigned int min_max_angle;
  unsigned int preserve_ridges;
  float flat_rate, dense_rate;
  float max_max_normal_angle;
  ActionItem* next;
};


// Constructor of ActionItem
ActionItem* new_ActionItem()
{
  ActionItem* action_item;
  action_item = (ActionItem*)malloc(sizeof(ActionItem));
  
  // Assign default values
  action_item->action               = no_action;
  action_item->max_iter             = MAX_ITER;
  action_item->num_iter             = NUM_ITER;
  action_item->max_min_angle        = MAX_MIN_ANGLE;
  action_item->min_max_angle        = MIN_MAX_ANGLE;
  action_item->preserve_ridges      = PRESERVE_RIDGES;
  action_item->dense_rate           = DENSE_RATE;
  action_item->flat_rate            = FLAT_RATE;
  action_item->max_max_normal_angle = MAX_MAX_NORMAL_ANGLE;
  action_item->next = NULL;
  return action_item;
}

// Destructor of ActionItem
void delete_ActionItem(ActionItem* action_list)
{
  ActionItem* tmp_list;
  while (action_list != NULL) {
    tmp_list = action_list->next;
    free(action_list);
    action_list = tmp_list;
  }
}

// Print usage of program
void usage()
{
  printf("");
  printf("Usage: ImproveSurfMesh [options...] file\n"
         "\n"
         " file - A surface mesh given in .off format\n"
         "\n"
	 "Several options can be combined. The order of the options will deside \n"
	 "in what order the mesh improvements are applied.\n"
	 "\n");
  printf("Options:\n");
  printf(" --smooth [MAX_ITER [PRESERVE_RIDGES [MAX_MIN_ANGLE MIN_MAX_ANGLE]]]\n"
         "      Improve the quality of the surface mesh using a angle-based\n"
         "      method. Geometrical features are preserved using a local structure\n"
         "      tensor. The algorithm is run until MAX_ITER is reached or the\n"
         "      angle conditions are saticefied.\n"
	 "\n"
         "      MAX_ITER[%d]: Maximal iteration of smooth cycles.\n"
	 "\n"
	 "      PRESERVE_RIDGES[%d]: Do not flip edges when it might destroy ridges.\n"\
	 "\n"
         "      MAX_MIN_ANGLE[%d]: The maximal value of the minimal angle.\n"
	 "\n"
	 "      MIN_MAX_ANGLE[%d]: The minimal value of the maximal angle.\n"
	 "\n", 
	 MAX_ITER, PRESERVE_RIDGES, MAX_MIN_ANGLE, MIN_MAX_ANGLE);
  printf(" --normal-smooth\n"
         "      Improve the quality of the surface mesh by smoothing the vertex\n"
         "      normals. Adviceable to do after a coarse - smooth cycle.\n"
	 "\n");
  printf(" --coarse-flat [RATE]\n"
         "      Remove vertices in flat areas\n"
	 "\n"
         "      RATE[%0.3f]: A cutoff value for deciding if a vertex will be removed.\n"
         "         A larger number will remove more vertices.\n" 
	 "\n", 
	 FLAT_RATE);
  printf(" --coarse-dense [NUM_ITER [RATE]]\n"
         "      Remove vertices in dense areas.\n"
	 "\n"
         "      NUM_ITER[%d]: Number of coarsening cycles.\n"
	 "\n"
         "      RATE[%0.1f]: A cutoff value for deciding if a vertex will be removed.\n"
         "         A larger number will remove more vertices.\n" 
	 "\n"
	 "\n", 
	 NUM_ITER, DENSE_RATE);
  printf(" --refine\n"
         "      Uniformly refine mesh\n"
	 "\n");
  printf(" --correct-normals\n"
         "      Correct the face normals so they all point outward.\n"
	 "\n");
  printf(" --active-sites FILE \n"
         "      A file with active sites which will not be coarsed.\n" 
	 "      See README for specifications.\n"
	 "\n"); 
  printf(" --stats [MAX_MIN_ANGLE MIN_MAX_ANGLE]\n"
         "      Print quality statistics of a surface mesh.\n"
	 "\n"
         "      MAX_MIN_ANGLE[%d]: The maximal value of the minimal angle.\n"
	 "\n"
	 "      MIN_MAX_ANGLE[%d]: The minimal value of the maximal angle.\n"
	 "\n", 
	 MAX_MIN_ANGLE, MIN_MAX_ANGLE);
  printf(" --help\n"
         "      Print this help.\n");
}

void option_error()
{
  usage();
  exit(1);
}

// Small helper function to increase the arg item variable
void increase_arg_item(int& arg_item, int argc) 
{
  arg_item++;
  if (arg_item == argc)
    option_error();
}

// Parse the command line arguments
void parse_comand_line(int argc, char *argv[], char*& input_filename,
		       char*& input_sites, ActionItem*& action_list)
{
  // Init local variables
  ActionItem* action_item = NULL;
  int arg_item = 1;
  
  // If not options given
  if (argc == arg_item)
    option_error();

  // Parse command line
  while (arg_item < argc){
    
    // Parse --active-sites option
    // ------------------------------------------------------------------
    if (strcmp (argv[arg_item], "--active-sites") == 0){
      increase_arg_item(arg_item, argc);

      // Get active site file name
      input_sites = argv[arg_item];
      increase_arg_item(arg_item, argc);
      
      continue;
    }
    
    // If the next argument is an action option
    if  (strncmp(argv[arg_item], "--", 2) == 0){
      
      // If it is the first action item created
      if (action_item == NULL){
	action_item = new_ActionItem();
	action_list = action_item;
      }
      else {
	action_item->next = new_ActionItem();
	action_item = action_item->next;
      }
    }
    
    // If the next argument is a letter assume input file name
    else if ( isalpha(argv[arg_item][0])){
      input_filename = argv[arg_item];
      arg_item++;
      if (arg_item != argc)
	option_error();
      break;
    }
    
    // Parse --help
    // ------------------------------------------------------------------
    if (strcmp (argv[arg_item], "--help") == 0)
      option_error();
    
    // Parse --refine option
    // ------------------------------------------------------------------
    if (strcmp (argv[arg_item], "--refine") == 0){
      action_item->action = refine;
      increase_arg_item(arg_item, argc);
      continue;
    }
    
    // Parse --normal-smooth option
    // ------------------------------------------------------------------
    if (strcmp (argv[arg_item], "--normal-smooth") == 0){
      action_item->action = normal_smooth;
      increase_arg_item(arg_item, argc);
      continue;
    }
    
    // Parse --correct-normals option
    // ------------------------------------------------------------------
    if (strcmp (argv[arg_item], "--correct-normals") == 0){
      action_item->action = correct_normals;
      increase_arg_item(arg_item, argc);
      continue;
    }
    
    // Parse --smooth option
    // ------------------------------------------------------------------
    if (strcmp (argv[arg_item], "--smooth") == 0){

      action_item->action = smooth;
      increase_arg_item(arg_item, argc);

      // If no values passed
      if  (strncmp(argv[arg_item], "--", 2) == 0 || isalpha(argv[arg_item][0]))
	continue;
      
      // Parse max iter argument
      action_item->max_iter = atoi(argv[arg_item]);
      increase_arg_item(arg_item, argc);
      
      // Check value
      if (action_item->max_iter < 1){
	printf("*** Argument error: Expected a positive value for the \'MAX_ITER\' argument.\n\n");
	exit(1);
      }

      // If no other values passed
      if (strncmp(argv[arg_item], "--", 2) == 0 || isalpha(argv[arg_item][0]))
	continue;
      
      // Parse flip edges argument
      action_item->preserve_ridges = atoi(argv[arg_item]);
      increase_arg_item(arg_item, argc);
      
      // Check value
      if (action_item->preserve_ridges != 1 && action_item->preserve_ridges != 0){
	printf("*** Argument error: Expected \'0\' or \'1\' for the \'PRESERVE_RIDGES\' argument.\n\n");
	exit(1);
      }

      // If no other values passed
      if (strncmp(argv[arg_item], "--", 2) == 0 || isalpha(argv[arg_item][0]))
	continue;
      
      // Parse max min angle
      action_item->max_min_angle = atoi(argv[arg_item]);
      increase_arg_item(arg_item, argc);
      
      // Check value
      if (action_item->max_min_angle < 1 || action_item->max_min_angle > 50){
	printf("*** Argument error: Expected a value between 1 and 50 for the \'MAX_MIN_ANGLE\' argument.\n\n");
	exit(1);
      }
      
      // Check that we get another value for MIN_MAX_ANGLE
      if (strncmp(argv[arg_item], "--", 2) == 0 || isalpha(argv[arg_item][0])){
	printf("*** Argument error: Expected a value for the \'MIN_MAX_ANGLE\' argument.\n\n");
	exit(1);
      }
      
      // Parse min max angle
      action_item->min_max_angle = atoi(argv[arg_item]);
      increase_arg_item(arg_item, argc);
      
      // Check value
      if (action_item->min_max_angle < 70 || action_item->min_max_angle > 179){
	printf("*** Argument error: Expected a value between 70 and 179 for the \'MIN_MAX_ANGLE\' argument.\n\n");
	exit(1);
      }
      
      continue;
    }

    // Parse --coarse-dense option
    // ------------------------------------------------------------------
    if (strcmp (argv[arg_item], "--coarse-dense") == 0){
      action_item->action = coarse_dense;
      increase_arg_item(arg_item, argc);
      
      // If no values passed
      if  (strncmp(argv[arg_item], "--", 2) == 0 || isalpha(argv[arg_item][0]))
	continue;
      
      // Parse num_iter argument
      action_item->num_iter = atoi(argv[arg_item]);
      increase_arg_item(arg_item, argc);
      
      // Check value
      if (action_item->num_iter <= 0){
	printf("*** Argument error: Expected a positive integer for the \'NUM_ITER\' argument.\n\n");
	exit(1);
      }
      
      // Check if we are finnished parsing coarsning parameters
      if (strncmp(argv[arg_item], "--", 2) == 0 || isalpha(argv[arg_item][0]))
	continue;
      
      // Parse dense_rate
      action_item->dense_rate = atof(argv[arg_item]);
      increase_arg_item(arg_item, argc);
      
      // Check value
      if (action_item->dense_rate < 0){
	printf("*** Argument error: Expected a positive value for the \'RATE\' argument.\n\n");
	exit(1);
      }
	
      continue;
    }
    
    // Parse --coarse-flat option
    // ------------------------------------------------------------------
    if (strcmp (argv[arg_item], "--coarse-flat") == 0){
      action_item->action = coarse_flat;
      increase_arg_item(arg_item, argc);
      
      // If no values passed
      if  (strncmp(argv[arg_item], "--", 2) == 0 || isalpha(argv[arg_item][0]))
	continue;
      
      // Parse flat_rate
      action_item->flat_rate = atof(argv[arg_item]);
      increase_arg_item(arg_item, argc);
      
      // Check value
      if (action_item->flat_rate < 0){
	printf("*** Argument error: Expected a positive value for the \'RATE\' argument.\n\n");
	exit(1);
      }
	
      continue;
    }
    
    // Parse --stats option
    // ------------------------------------------------------------------
    if (strcmp (argv[arg_item], "--stats") == 0){

      action_item->action = stats;
      increase_arg_item(arg_item, argc);

      // If no values passed
      if  (strncmp(argv[arg_item], "--", 2) == 0 || isalpha(argv[arg_item][0]))
	continue;
      
      // Parse max min angle
      action_item->max_min_angle = atoi(argv[arg_item]);
      increase_arg_item(arg_item, argc);
      
      // Check value
      if (action_item->max_min_angle < 1 || action_item->max_min_angle > 50){
	printf("*** Argument error: Expected a value between 1 and 50 for the \'MAX_MIN_ANGLE\' argument.\n\n");
	exit(1);
      }
      
      // Check that we get another value for MIN_MAX_ANGLE
      if (strncmp(argv[arg_item], "--", 2) == 0 || isalpha(argv[arg_item][0])){
	printf("*** Argument error: Expected a value for the \'MIN_MAX_ANGLE\' argument.\n\n");
	exit(1);
      }
      
      // Parse min max angle
      action_item->min_max_angle = atoi(argv[arg_item]);
      increase_arg_item(arg_item, argc);
      
      // Check value
      if (action_item->min_max_angle < 70 || action_item->min_max_angle > 179){
	printf("*** Argument error: Expected a value between 70 and 179 for the \'MIN_MAX_ANGLE\' argument.\n\n");
	exit(1);
      }
      
      continue;
    }

    printf("Unknown argument: \'%s\'\n", argv[arg_item]);
    exit(1);

  }

  // Check for no action options
  if (action_item == NULL) {
    printf("No action option choosen. Nothing to be done.\n");
    exit(1);
  }
  
}

/*
 * ***************************************************************************
 * Routine:  ReadSurfMeshFile
 *
 * Author:   Johan Hake (hake.dev@gmail.com), Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Read a surface mesh file stored in OFF format
 * ***************************************************************************
 */
SurfaceMesh* ReadSurfMeshFile(char* input_name, bool read_suffix_count, 
			      char* basename, int& suffix_count)
{
  char filename[256];

  char* basename_end_pointer;
  char* suffix_end_pointer;
  unsigned int suffix_index;
  unsigned int basename_index;
  int suffix_offset;
  SurfaceMesh* surfmesh;

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

  // Check if input has _suffix_XX in its name
  if (read_suffix_count){
    suffix_end_pointer = strstr(filename, OUT_SUFFIX);
    if (suffix_end_pointer != NULL){
      // If so update basename and get suffix_count
      suffix_index = suffix_end_pointer-filename;
      strncpy (basename, filename, suffix_index);
      basename[suffix_index] = '\0';
      suffix_offset = (basename_index - (suffix_index + OUT_SUFFIX_LENGTH));
    
      if (suffix_offset == 1 || suffix_offset == 2)
	sscanf(basename_end_pointer-suffix_offset, "%d", &suffix_count);
    }
  }

  // Load surface meshes in OFF format
  return SurfaceMesh_readOFF(input_name);
}

/*
 * ***************************************************************************
 * Routine:  ImproveSurfaceMesh
 *
 * Author:   Johan Hake (hake.dev@gmail.com), Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Improve a surface mesh
 * ***************************************************************************
 */
int ImproveSurfaceMesh(char* input_name, char* input_site, ActionItem* action_list)
{
  FILE *fout;
  time_t t0, t1, t2;
  SurfaceMesh* surfmesh = NULL;
  char filename[256], basename[256];
  unsigned int num_spheres = 0;
  ATOM *sphere_list = NULL;
  unsigned int* sphere_markers = NULL;
  //float* area_constraint_list = NULL; // Not used
  int suffix_count = -1;
  ActionItem* action = action_list;
  float acos_max_max_normal_angle;
  unsigned int i, j;

  (void)time(&t0);
  
  // Read any active sites from file
  if (input_site != NULL)
    ReadActiveSiteFile(input_site, num_spheres,
  //                     sphere_list, sphere_markers, area_constraint_list);
                       sphere_list, sphere_markers);
  // Read surface mesh file together with basename and suffix count
  surfmesh = ReadSurfMeshFile(input_name, true, basename, suffix_count);

  // Create neighborhood information 
  SurfaceMesh_createNeighborlist(surfmesh);

  // Main action loop
  while (action != NULL){

    switch(action->action){
      
      // Refine mesh
      // ------------------------------------------------------------------
      case refine: {

	(void)time(&t1);
	
	printf("\n\n");
	printf("Refine mesh:\n");

	SurfaceMesh_refine(surfmesh);

	printf("After refinement:              \n"
	       "   vertices: %d --- simplices: %d\n", surfmesh->num_vertices, surfmesh->num_faces);
      
	(void)time(&t2);
	printf("------------------------------------------\n");  
	printf("Time to refine surface mesh: %4.1f seconds. \n", fabs(t2-t1));  
	
	// Set changed flag
	mesh_changed = true;

	break;
      }
	
      // Smooth mesh
      // ------------------------------------------------------------------
      case smooth: {
	(void)time(&t1);
	printf("\n\n");
	printf("Smooth mesh:\n");
	SurfaceMesh_smooth(surfmesh, action->max_min_angle, 
			   action->min_max_angle, action->max_iter, 
			   action->preserve_ridges == 1);
	(void)time(&t2);
	printf("------------------------------------------\n");  
	printf("Time to smooth surface mesh: %4.1f seconds. \n", fabs(t2-t1));  

	// Set changed flag
	mesh_changed = true;

	break;
      }
    
      // Normal smooth mesh
      // ------------------------------------------------------------------
      case normal_smooth: {
	(void)time(&t1);
	printf("\n\n");
	printf("Normal smooth mesh:\n");
	SurfaceMesh_normalSmooth(surfmesh);
	(void)time(&t2);
	printf("-------------------------------------------------\n");  
	printf("Time to normal smooth surface mesh: %4.1f seconds.\n", fabs(t2-t1));  

	// Set changed flag
	mesh_changed = true;

	break;
      }

      // Correct normals of surface mesh
      // ------------------------------------------------------------------
      case correct_normals: {
	(void)time(&t1);
	printf("\n\n");
	printf("Correct normals:\n");
	SurfaceMesh_correctNormals(surfmesh);
	(void)time(&t2);
	printf("---------------------------------------------------\n");  
	printf("Time to correct surface mesh normals: %4.1f seconds.\n", fabs(t2-t1));  

	// Set changed flag
	mesh_changed = true;

	break;
      }

      // Coarse mesh
      // ------------------------------------------------------------------
      case coarse_dense: {
	(void)time(&t1);

	printf("\n\n");
	printf("Coarse dense part of mesh:\n");

	// Do the actuall coarsening
	for (i=0; i<action->num_iter; i++){
	  
	  // Assign active sites
	  SurfaceMesh_assignActiveSites(surfmesh, sphere_list, num_spheres, 
					sphere_markers);
	  
	  // Do the actuall coarsening
	  SurfaceMesh_coarse(surfmesh, action->dense_rate, 0, 10, -1);
	  if (action->num_iter==1)
	    printf("After coarsening:        \n   vertices: %d --- simplices: %d\n", surfmesh->num_vertices, surfmesh->num_faces);
	  else
	    printf("%d After coarsening:        \n   vertices: %d --- simplices: %d\n", i, surfmesh->num_vertices, surfmesh->num_faces);
	}
	
	(void)time(&t2);
	printf("------------------------------------------\n");  
	printf("Time to coarse surface mesh: %4.1f seconds.\n", fabs(t2-t1));  

	// Set changed flag
	mesh_changed = true;

	break;
      }

      // Coarse mesh
      // ------------------------------------------------------------------
      case coarse_flat: {
	(void)time(&t1);

	printf("\n\n");
	printf("Coarse flat part of mesh:\n");

	// Assign active sites
	SurfaceMesh_assignActiveSites(surfmesh, sphere_list, num_spheres,
				      sphere_markers);
	
	// Do the actuall coarsening
	SurfaceMesh_coarse(surfmesh, action->flat_rate, 1, 0, -1);
	
	printf("After coarsening:        \n   vertices: %d --- simplices: %d\n", surfmesh->num_vertices, surfmesh->num_faces);
      
	(void)time(&t2);
	printf("------------------------------------------\n");  
	printf("Time to coarse surface mesh: %4.1f seconds.\n", fabs(t2-t1));  

	// Set changed flag
	mesh_changed = true;

	break;
      }

      // Get stats
      // ------------------------------------------------------------------
      case stats: {
	float min_angle, max_angle;
	int num_small, num_large;

	printf("\n\n");
	printf("Surface mesh statistics:\n");
	SurfaceMesh_getMinMaxAngles(surfmesh ,&min_angle, &max_angle, &num_small,
				    &num_large, action->max_min_angle, 
				    action->min_max_angle);

	printf("    min_angle: %f - max_angle: %f - "
	       "smaller-than-%d: %d - larger-than-%d: %d\n", 
	       min_angle, max_angle, action->max_min_angle, num_small, 
	       action->min_max_angle, num_large);

	break;
      }
    
    }
    
    // The next action
    action = action->next;

  }

  // Ouput mesh if it has changed
  if (mesh_changed)
  {

    // Output surface mesh
    sprintf(filename, "%s%s%d.off", basename, OUT_SUFFIX, ++suffix_count);

    printf("\n\n");
    printf("Writing results to: %s%s%d.off\n", basename, OUT_SUFFIX, suffix_count);

    // Write surface mesh to file
    SurfaceMesh_writeOFF(surfmesh, filename);
  
    (void)time(&t2);
    
    printf("\n");
    printf("----------------------------------------\n");  
    printf("Total time improving mesh: %4.1f seconds. \n\n", fabs(t2-t0));  
  }

  // Release memory
  SurfaceMesh_dtor(surfmesh);
  
  // If any active sites have been read
  if (num_spheres){
    free(sphere_list);
    free(sphere_markers);
  }
  
  return(0);
}

// The main function
// ------------------------------------------------------------------
int main(int argc, char *argv[])
{
  char* input_sites = NULL;
  char* input_filename = NULL;
  ActionItem* action_list;

  // Parse command line
  parse_comand_line(argc, argv, input_filename, input_sites, action_list);

  // Call ImproveSurfaceMesh routine
  ImproveSurfaceMesh(input_filename, input_sites, action_list);

  // Clean action list
  delete_ActionItem(action_list);

  return(0);
}
