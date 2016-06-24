/*
 * ***************************************************************************
 * GAMER = < Geometry-preserving Adaptive MeshER >
 * Copyright (C) 2007-2010 -- Michael Holst and Johan Hake
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
 * File:     Triangulate.C    < ... >
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Triangulate (using Triangle) a connected list of vertices
 * ****************************************************************************
 */

#include <gamer/biom.h>
#include <gamer/triangle.h>
#include "gamercf.h"

triangulateio* init_triangulateio()
{
  triangulateio* io = new triangulateio;
  io->pointlist  = (double *) NULL;
  io->pointattributelist = (double *) NULL;
  io->pointmarkerlist = (int *) NULL;
  io->numberofpoints = 0;
  io->numberofpointattributes = 0;

  io->trianglelist = (int *) NULL;
  io->triangleattributelist = (double *) NULL;
  io->trianglearealist = (double *) NULL;
  io->neighborlist = (int *) NULL;
  io->numberoftriangles = 0;
  io->numberofcorners = 3;
  io->numberoftriangleattributes = 0;

  io->segmentlist = (int *) NULL;
  io->segmentmarkerlist = (int *) NULL;
  io->numberofsegments = 0;

  io->holelist = (double *) NULL;
  io->numberofholes = 0;

  io->regionlist = (double *) NULL;
  io->numberofregions = 0;      

  io->edgelist = (int *) NULL;
  io->edgemarkerlist = (int *) NULL;
  io->normlist = (double *) NULL;
  io->numberofedges = 0;
}

SurfaceMesh* SurfaceMesh_triangulate(REAL* pointlist, int numberofcoordinates,
				     char* triangle_params)
{
  // Start with empty triangle structures
  triangulateio* in = init_triangulateio();
  triangulateio* out = init_triangulateio();

  // Attach the given points
  in->pointlist = pointlist;
  in->numberofpoints = numberofcoordinates/2;
  
  // Create segments list
  in->numberofsegments = numberofcoordinates/2;
  in->segmentlist = new int[numberofcoordinates];

  // Fill segment list
  for (int i = 0; i < numberofcoordinates/2; i++)
  {
    // First vertices is always i
    in->segmentlist[i*2] = i;
    
    // Special case for last segment
    if (i == numberofcoordinates/2-1)
      in->segmentlist[i*2+1] = 0;
    else
      in->segmentlist[i*2+1] = i+1;
  }
  
  // Call triangulate from triangle
  triangulate(triangle_params, in, out, (struct triangulateio *) NULL);
  
  // Initiate an empty SurfaceMesh
  SurfaceMesh* surfmesh = SurfaceMesh_ctor(out->numberofpoints, out->numberoftriangles);
  
  // Grab vertex coordinates from the triangulated mesh
  for (int n=0; n < out->numberofpoints; n++)
  {
    surfmesh->vertex[n].x = out->pointlist[n*2];
    surfmesh->vertex[n].y = out->pointlist[n*2+1];
    surfmesh->vertex[n].z = 0.0;
  }

  // Grab connectivity from the triangulated mesh
  for (int n=0; n < out->numberoftriangles; n++)
  {
    surfmesh->face[n].a = out->trianglelist[n*3];
    surfmesh->face[n].b = out->trianglelist[n*3+1];
    surfmesh->face[n].c = out->trianglelist[n*3+2];
  }
  
  // Clean up
  free(out->pointlist);
  free(out->trianglelist);
  delete[] in->segmentlist;
  delete out;
  delete in;
  return surfmesh;
}
