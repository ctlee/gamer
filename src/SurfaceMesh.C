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
 * File:     SurfaceMeshOld.C    < ... >
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Create and destroy SurfaceMeshOld data
 * ****************************************************************************
 */

#include "biom.h"
#include "SurfaceMeshOld.h"
#include <vector>
#include <algorithm>
#include <utility>
#include <cmath>

// Declare internal GAMer methods
EIGENVECT GetEigenVector(SurfaceMeshOld *,
                         int,
                         FLTVECT *,
                         float *);

/*
 * ***************************************************************************
 * Routine:  SurfaceMeshOld
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 * c++ Revisions: John Moody (brogan@gmail.com)
 *
 * Purpose:  Create a surface mesh instance
 * ***************************************************************************
 */
SurfaceMeshOld::SurfaceMeshOld(unsigned int num_vertices, unsigned int num_faces)
{
    this->num_vertices = num_vertices;
    this->num_faces    = num_faces;

    if (num_vertices)
    {
        vertex = (FLTVECT *)malloc(sizeof(FLTVECT) * num_vertices);
    }
    else
    {
        vertex = NULL;
    }

    if (num_faces)
    {
        face = (INT3VECT *)malloc(sizeof(INT3VECT) * num_faces);
    }
    else
    {
        face = NULL;
    }

    neighbor      = NULL;
    neighbor_list = NULL;
    avglen        = 0.;
    min[0]        = 0.; min[1] = 0.; min[2] = 0.;
    max[0]        = 0.; max[1] = 0.; max[2] = 0.;

    // Initialize SurfaceMeshOld structures
    for (int n = 0; n < num_vertices; n++)
    {
        FLTVECT& vert = vertex[n];
        vert.x   = vert.y = vert.z = vert.m = 0;
        vert.sel = true;
    }

    for (int n = 0; n < num_faces; n++)
    {
        INT3VECT& f = face[n];
        f.a   = f.b = f.c = f.m = 0;
        f.sel = true;
    }

    // Initialize domain data
    closed                = true;
    _marker                = 1;
    volume_constraint     = 100;
    use_volume_constraint = false;
    hole                  = false;
}

/*
 * ***************************************************************************
 * Routine:  ~SurfaceMeshOld()
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 * c++ Revisions: John Moody (brogan@gmail.com)
 *
 * Purpose:  Release memory from a SurfaceMeshOld instance
 * ***************************************************************************
 */
SurfaceMeshOld::~SurfaceMeshOld()
{
    // Free allocated memory
    if (vertex)
    {
        free(vertex);
        vertex = 0;
    }

    if (face)
    {
        free(face);
        face = 0;
    }

    // Destroy neighbor_list
    destroyNeighborlist();
}

/*
 * ***************************************************************************
 * SubRoutine:  SurfaceMeshOld::createNeighborList
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com) and Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Create a neighbor list from a surface mesh
 * ***************************************************************************
 */
void SurfaceMeshOld::createNeighborlist()
{
    int n, a0, b0;
    int a, b, c, d;
    NPNT3 *first_ngr, *second_ngr, *tmp_ngr;
    bool   closed;

    if (neighbor_list)
    {
        return;
    }

    // Destroy any exsisting neighborlist
    destroyNeighborlist();

    // Create an array of NPNT3, used to store
    neighbor_list = (NPNT3 **)malloc(sizeof(NPNT3 *) * num_vertices);

    // Initialize the neighbor list
    for (n = 0; n < num_vertices; n++)
    {
        neighbor_list[n] = NULL;

        // Default mark all vertices for deletion
        vertex[n].m = -1;
    }

    // Iterate over the faces and collect line segments (a, b) and its connection
    // to a face (c). Save the line segment so it forms a counter clockwise triangle
    // with the origin vertex
    int num_connected = 0;

    for (n = 0; n < num_faces; n++)
    {
        a = face[n].a;
        b = face[n].b;
        c = face[n].c;

        if ((a == b) || (b == c) || (a == c))
        {
            printf("Face %d include vertices with same indices (%d, %d, %d).\n", n, a, b, c);
        }

        first_ngr        = (NPNT3 *)malloc(sizeof(NPNT3));
        first_ngr->a     = b;
        first_ngr->b     = c;
        first_ngr->c     = n;
        first_ngr->next  = neighbor_list[a];
        neighbor_list[a] = first_ngr;

        // Mark vertex as connected
        if (vertex[a].m < 0)
        {
            vertex[a].m    = 0;
            num_connected += 1;
        }

        first_ngr        = (NPNT3 *)malloc(sizeof(NPNT3));
        first_ngr->a     = c;
        first_ngr->b     = a;
        first_ngr->c     = n;
        first_ngr->next  = neighbor_list[b];
        neighbor_list[b] = first_ngr;

        // Mark vertex as connected
        if (vertex[b].m < 0)
        {
            vertex[b].m    = 0;
            num_connected += 1;
        }

        first_ngr        = (NPNT3 *)malloc(sizeof(NPNT3));
        first_ngr->a     = a;
        first_ngr->b     = b;
        first_ngr->c     = n;
        first_ngr->next  = neighbor_list[c];
        neighbor_list[c] = first_ngr;

        // Mark vertex as connected
        if (vertex[c].m < 0)
        {
            vertex[c].m    = 0;
            num_connected += 1;
        }
    }

    // Check if there are vertices which are not connect to any face
    if (num_connected < num_vertices)
    {
        destroyNeighborlist();

        // Remove unconnected vertices
        removeUnconnectedVertices();

        // Re-create neighbors
        createNeighborlist();

        return;
    }

    // Order the neighbors so they are connected counter clockwise
    for (n = 0; n < num_vertices; n++)
    {
        a0        = b0 = -1;
        first_ngr = neighbor_list[n];
        c         = first_ngr->a;
        d         = first_ngr->b;

        while (first_ngr != NULL)
        {
            a = first_ngr->a;
            b = first_ngr->b;

            second_ngr = first_ngr->next;

            while (second_ngr != NULL)
            {
                a0 = second_ngr->a;
                b0 = second_ngr->b;

                if ((a0 == b) && (b0 != a))
                {
                    tmp_ngr = first_ngr;

                    while (tmp_ngr != NULL)
                    {
                        if (tmp_ngr->next == second_ngr)
                        {
                            tmp_ngr->next = second_ngr->next;
                            break;
                        }
                        tmp_ngr = tmp_ngr->next;
                    }
                    tmp_ngr          = first_ngr->next;
                    first_ngr->next  = second_ngr;
                    second_ngr->next = tmp_ngr;
                    break;
                }

                second_ngr = second_ngr->next;
            }

            first_ngr = first_ngr->next;
        }

        // Check that the neighbor list is connected
        tmp_ngr = neighbor_list[n];

        closed = true;

        while (tmp_ngr->next != NULL)
        {
            // Check that we are connected
            if (tmp_ngr->b != tmp_ngr->next->a)
            {
                if (closed)
                {
                    printf("Polygons connected to vertex %d are not closed (interupted):"
                           " (%.2f, %.2f, %.2f)\n", n, vertex[n].x,
                           vertex[n].y, vertex[n].z);
                }

                // Do not bail, just register the vertex to not be done anything with
                vertex[n].sel = false;

                closed = false;
            }

            // Step one face forward
            tmp_ngr = tmp_ngr->next;
        }

        // Check if the list forms a closed ring
        if (closed && (b0 != c))
        {
            printf("Polygons connected to vertex %d are not closed (not closed):"
                   " (%.2f, %.2f, %.2f)\n", n, vertex[n].x,
                   vertex[n].y, vertex[n].z);

            // Do not bail, just register the vertex to not be done anything with
            vertex[n].sel = false;

            closed = false;
        }

        if (!closed)
        {
            closed = false;
        }
    }
}

/*
 * ***************************************************************************
 * SubRoutine:  SurfaceMeshOld::destroyNeighborlist
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 * c++ revisions: John Moody (brogan@gmail.com)
 *
 * Purpose:  Release memory from a neighborlist
 * ***************************************************************************
 */
void SurfaceMeshOld::destroyNeighborlist()
{
    unsigned int n;
    NPNT3 *first_ngr = NULL;
    NPNT3 *tmp_ngr   = NULL;

    if (neighbor_list != NULL)
    {
        // Release the single neighbors
        for (n = 0; n < num_vertices; n++)
        {
            first_ngr = neighbor_list[n];

            while (first_ngr != NULL)
            {
                tmp_ngr = first_ngr->next;
                free(first_ngr);
                first_ngr = tmp_ngr;
            }
        }

        // Free the array of pointers
        free(neighbor_list);
        neighbor_list = NULL;
    }
}

/*
 * ***************************************************************************
 * SubRoutine:  SurfaceMeshOld::flipNormals
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Release memory from a neighborlist
 * ***************************************************************************
 */
void SurfaceMeshOld::flipNormals()
{
    // First destroy neighborlist
    destroyNeighborlist();

    // Iterate over faces and flip the vertices
    for (int n = 0; n < num_faces; n++)
    {
        const int tmp = face[n].a;
        face[n].a = face[n].b;
        face[n].b = tmp;
    }
}

/*
 * ***************************************************************************
 * SubRoutine: SurfaceMeshOld::remove_unconnected_vertices
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Remove vertices which are not connected to any face
 * ***************************************************************************
 */
void SurfaceMeshOld::removeUnconnectedVertices()
{
    // Collect statistics
    int num_to_be_removed = 0;

    std::vector<int> to_be_removed(num_vertices);

    for (int n = 0; n < num_vertices; n++)
    {
        if (vertex[n].m < 0)
        {
            num_to_be_removed++;
        }
        to_be_removed[n] = num_to_be_removed;
    }

    printf("Removing %d vertices.\n", num_to_be_removed);

    // Move vertices forward
    for (int n = 0; n < num_vertices; n++)
    {
        // If a vertex is to be removed
        if (((n == 0) && (to_be_removed[n] != 0)) ||
            ((n != 0) && (to_be_removed[n - 1] != to_be_removed[n])))
        {
            continue;
        }

        // Move vertices forward
        vertex[n - to_be_removed[n]].x   = vertex[n].x;
        vertex[n - to_be_removed[n]].y   = vertex[n].y;
        vertex[n - to_be_removed[n]].z   = vertex[n].z;
        vertex[n - to_be_removed[n]].sel = vertex[n].sel;
        vertex[n - to_be_removed[n]].m   = vertex[n].m;
    }

    // Fix face offset
    for (int n = 0; n < num_faces; n++)
    {
        face[n].a = face[n].a - to_be_removed[face[n].a];
        face[n].b = face[n].b - to_be_removed[face[n].b];
        face[n].c = face[n].c - to_be_removed[face[n].c];
    }

    // Adjust num_vertices
    num_vertices -= num_to_be_removed;
}

/*
 * ***************************************************************************
 * SubRoutine: SurfaceMeshOld::deleteVertices
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Remove vertices which are marked with -1 as a vertex marker
 * ***************************************************************************
 */
void SurfaceMeshOld::deleteVertices()
{
    // Mark faces connected to vertices for deletion
    for (int n = 0; n < num_faces; n++)
    {
        if ((vertex[face[n].a].m < 0) ||
            (vertex[face[n].b].m < 0) ||
            (vertex[face[n].c].m < 0))
        {
            face[n].m = -1;
        }
    }

    // Delete marked faces
    deleteFaces();
}

/*
 * ***************************************************************************
 * SubRoutine: SurfaceMeshOld::deleteFaces
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Remove faces which are marked with -1 as a face marker
 * ***************************************************************************
 */
void SurfaceMeshOld::deleteFaces()
{
    // Iterate over vertices and mark all for deletion
    for (int n = 0; n < num_vertices; n++)
    {
        vertex[n].m = -1;
    }

    // Delete faces connected to vertices
    int num_removed = 0;

    for (int n = 0; n < num_faces; n++)
    {
        // Check for removal of face
        if (face[n].m < 0)
        {
            num_removed += 1;
        }
        else
        {
            // If any previous face has been marked for deletion
            if (num_removed > 0)
            {
                // Copy the face to a previous face
                face[n - num_removed].a   = face[n].a;
                face[n - num_removed].b   = face[n].b;
                face[n - num_removed].c   = face[n].c;
                face[n - num_removed].m   = face[n].m;
                face[n - num_removed].sel = face[n].sel;
            }

            // Un mark vertex for deletion
            vertex[face[n].a].m = 0;
            vertex[face[n].b].m = 0;
            vertex[face[n].c].m = 0;
        }
    }

    // Update the number of faces
    num_faces -= num_removed;
    removeUnconnectedVertices();
}

/*
 * ***************************************************************************
 * SubRoutine: SurfaceMeshOld::getCenterRadius
 *
 * Author:   Johan Hake (hake.dev@gmail.com)
 *
 * Purpose:  Return the center and radius of SurfaceMeshOld
 * ***************************************************************************
 */
ATOM SurfaceMeshOld::getCenterRadius()
{
    ATOM data;
    int  i;

    // Get the center and radius of the molecular surface mesh
    float mol_center_x = 0;
    float mol_center_y = 0;
    float mol_center_z = 0;
    float distance;

    for (i = 0; i < num_vertices; i++)
    {
        mol_center_x += vertex[i].x;
        mol_center_y += vertex[i].y;
        mol_center_z += vertex[i].z;
    }

    if (num_vertices > 0)
    {
        data.x = (float)(mol_center_x / (double)num_vertices);
        data.y = (float)(mol_center_y / (double)num_vertices);
        data.z = (float)(mol_center_z / (double)num_vertices);
    }
    else
    {
        printf("no data found ...\n");
        data.radius = 0.0;
        data.x      = 0.0;
        data.y      = 0.0;
        data.z      = 0.0;
        return data;
    }

    data.radius = 0;

    for (i = 0; i < num_vertices; i++)
    {
        distance = sqrt((vertex[i].x - data.x) * (vertex[i].x - data.x) +
                        (vertex[i].y - data.y) * (vertex[i].y - data.y) +
                        (vertex[i].z - data.z) * (vertex[i].z - data.z));

        if (distance > data.radius)
        {
            data.radius = distance;
        }
    }

    return data;
}

void SurfaceMeshOld::translate(float dx, float dy, float dz)
{
    int i;

    for (i = 0; i < num_vertices; i++)
    {
        vertex[i].x += dx;
        vertex[i].y += dy;
        vertex[i].z += dz;
    }
}

void SurfaceMeshOld::scale(float scale_x,
                        float scale_y, float scale_z)
{
    int i;

    for (i = 0; i < num_vertices; i++)
    {
        vertex[i].x *= scale_x;
        vertex[i].y *= scale_y;
        vertex[i].z *= scale_z;
    }
}

void SurfaceMeshOld::scale(float s)
{
    scale(s, s, s);
}

void SurfaceMeshOld::centeralize()
{
    ATOM center = getCenterRadius();

    translate(-center.x, -center.y, -center.z);
}

void SurfaceMeshOld::eigenvalues()
{
    EIGENVECT eigen_vect;
    FLTVECT   eigen_value;
    float     max_angle;
    int n;

    // Check if neighborlist is created
    if (!neighbor_list)
    {
        createNeighborlist();
    }

    if (neighbor_list == NULL)
    {
        printf("Could not create neighbor list some polygons might not be closed.\n");
        printf("Bailing out...\n");
        return;
    }

    for (n = 0; n < num_vertices; n++)
    {
        // If we have a vertex wich is not selected we continue
        if (!vertex[n].sel)
        {
            continue;
        }

        eigen_vect = GetEigenVector(this, n, &eigen_value, &max_angle);

        // FIXME: Fix a way to get this information out in some way...
        printf("%d: (%7.2f, %7.2f, %7.2f), %5.2f, %5.2f, %5.2f\n", n, vertex[n].x,
               vertex[n].y, vertex[n].z, eigen_value.x, eigen_value.y, eigen_value.z);
    }
}


void SurfaceMeshOld::splitMultipleConnectedSurfaces()
{
    int m, n, b0 = 0;
    int a, c = 0;
    NPNT3 *tmp_ngr, *last_ngr;
    bool   closed;

    // First stl library class... Let there be more...
    std::vector<int> non_manifold_vertices(0);

    // Re-create the neigborlist
    createNeighborlist();

    // Nothing to do
    if (this->isClosed())
    {
        return;
    }

    // Find non_manifold_vertices
    for (n = 0; n < num_vertices; n++)
    {
        // Check that the neighbor list is connected
        tmp_ngr = neighbor_list[n];

        closed = true;

        while (tmp_ngr->next != NULL)
        {
            // Check that we are connected
            if (tmp_ngr->b != tmp_ngr->next->a)
            {
                if (std::find(non_manifold_vertices.begin(), non_manifold_vertices.end(), n) ==
                    non_manifold_vertices.end())
                {
                    non_manifold_vertices.push_back(n);
                }

                // Do not bail, just register the vertex to not be done anything with
                vertex[n].sel = false;

                closed = false;
            }

            // Step one face forward
            tmp_ngr = tmp_ngr->next;
        }

        // Check if the list forms a closed ring
        if (closed && (b0 != c))
        {
            closed = false;

            // FIXME: Should we add non_manifold vertex for this case too?
            // FIXME: This code does not work at all. -JBM
        }
    }

    int new_number_of_vertices = num_vertices;

    printf("Experimental method for splitting multiple connected surfaces mesh.\n");
    printf("Old number of vertices: %d\n", num_vertices);

    // Scratch data
    std::vector<NPNT3 *> tmp_neighborgs(0);
    std::vector<FLTVECT> tmp_vertices(0);
    std::vector<int> not_fixed_verts(0);

    printf("Trying to repair miss connected mesh.\n");
    printf("Number of non manifold vertices: %d,\n", (int)non_manifold_vertices.size());

    for (int j = 0; j < non_manifold_vertices.size(); j++)
    {
        // Get faulty vertex
        n       = non_manifold_vertices[j];
        tmp_ngr = neighbor_list[n];
        a       = tmp_ngr->a;
        bool fixed      = true;
        int  num_first  = 0;
        int  num_second = 0;

        // Gather data for averaging vert coordinates of connected vertices
        FLTVECT first_vert, second_vert;
        first_vert.x  = first_vert.y = first_vert.z = 0;
        second_vert.x = second_vert.y = second_vert.z = 0;

        while (tmp_ngr->next != NULL)
        {
            num_first++;

            // Gather geometry information
            first_vert.x += vertex[tmp_ngr->a].x;
            first_vert.y += vertex[tmp_ngr->a].y;
            first_vert.z += vertex[tmp_ngr->a].z;

            // If we have closed a ring but are not finished
            if ((tmp_ngr->b != tmp_ngr->next->a) && (tmp_ngr->b == a))
            {
                // printf("Got one closed ring of %d faces for vert %d.\n", num_first, n);

                // Start a new neighbor ring
                last_ngr = tmp_ngr;
                tmp_ngr  = tmp_ngr->next;

                // Push back the start of the next ring
                tmp_neighborgs.push_back(tmp_ngr);

                c = tmp_ngr->a;

                // Walk the next ring of connected vertices
                while (tmp_ngr->next != NULL)
                {
                    num_second++;

                    // Gather geometry information
                    second_vert.x += vertex[tmp_ngr->a].x;
                    second_vert.y += vertex[tmp_ngr->a].y;
                    second_vert.z += vertex[tmp_ngr->a].z;

                    // The ring has a hole
                    if (tmp_ngr->b != tmp_ngr->next->a)
                    {
                        // If not already registered
                        // if (std::find(not_fixed_verts.begin(), not_fixed_verts.end(), n) ==
                        //	  not_fixed_verts.end())
                        //	printf("Could not fix vertex %d.\n", n);
                        // else
                        not_fixed_verts.push_back(n);

                        fixed             = false;
                        tmp_neighborgs[j] = NULL;
                        break;
                    }

                    tmp_ngr = tmp_ngr->next;
                }

                // Add last neighbor
                num_second++;

                // Gather geometry information
                second_vert.x += vertex[tmp_ngr->a].x;
                second_vert.y += vertex[tmp_ngr->a].y;
                second_vert.z += vertex[tmp_ngr->a].z;

                // The ring is not closed
                if (tmp_ngr->b != c)
                {
                    // If not already registered
                    if (std::find(not_fixed_verts.begin(), not_fixed_verts.end(), n) ==
                        not_fixed_verts.end())
                    {
                        printf("Could not fix vertex %d.\n", n);
                    }
                    else
                    {
                        not_fixed_verts.push_back(n);
                    }
                    fixed = false;
                    printf("Could not fix vertex %d.\n", n);
                    tmp_neighborgs[j] = NULL;
                }

                // If we manage to fix it
                if (fixed)
                {
                    // printf("Got a second closed ring of %d faces for vert %d.\n", num_second, n);

                    // Bump the total number of vertices
                    new_number_of_vertices += 1;

                    // Select vertex
                    vertex[n].sel = true;

                    // Register the vertices for addition
                    tmp_vertices.push_back(vertex[n]);
                    FLTVECT& tmp_vertex = tmp_vertices[tmp_vertices.size() - 1];

                    // Find mid point of vertex rings
                    first_vert.x  /= num_first;
                    first_vert.y  /= num_first;
                    first_vert.z  /= num_first;
                    second_vert.x /= num_second;
                    second_vert.y /= num_second;
                    second_vert.z /= num_second;

                    // Find midpoint between original vertex and vertex ring and
                    // update the two vertices with that point
                    // printf("old coordinate: (%.2f, %.2f, %.2f)\n", vertex[n].x,
                    //	   vertex[n].y, vertex[n].z);

                    vertex[n].x  = (vertex[n].x + first_vert.x) / 2;
                    vertex[n].y  = (vertex[n].y + first_vert.y) / 2;
                    vertex[n].z  = (vertex[n].z + first_vert.z) / 2;
                    tmp_vertex.x = (tmp_vertex.x + second_vert.x) / 2;
                    tmp_vertex.y = (tmp_vertex.y + second_vert.y) / 2;
                    tmp_vertex.z = (tmp_vertex.z + second_vert.z) / 2;


                    // printf("New coordinate: (%.2f, %.2f, %.2f)\n", vertex[n].x,
                    //   vertex[n].y, vertex[n].z);
                    // printf("New coordinate2: (%.2f, %.2f, %.2f)\n", tmp_vertex.x,
                    //	   tmp_vertex.y, tmp_vertex.z);

                    // Update neighborlist of neighbors
                    tmp_ngr = tmp_neighborgs[j];

                    while (tmp_ngr != NULL)
                    {
                        // Update face information
                        if (face[tmp_ngr->c].a == n)
                        {
                            // printf("Changing connection for face %d.a (%d->%d)\n", tmp_ngr->c,
                            //       n, new_number_of_vertices-1);
                            face[tmp_ngr->c].a = new_number_of_vertices - 1;
                        }

                        if (face[tmp_ngr->c].b == n)
                        {
                            // printf("Changing connection for face %d.b (%d->%d)\n", tmp_ngr->c,
                            //       n, new_number_of_vertices-1);
                            face[tmp_ngr->c].b = new_number_of_vertices - 1;
                        }

                        if (face[tmp_ngr->c].c == n)
                        {
                            // printf("Changing connection for face %d.c (%d->%d)\n", tmp_ngr->c,
                            //       n, new_number_of_vertices-1);
                            face[tmp_ngr->c].c = new_number_of_vertices - 1;
                        }

                        tmp_ngr = tmp_ngr->next;
                    }

                    break;
                }
            }

            // Step one face forward
            tmp_ngr = tmp_ngr->next;

            if (tmp_ngr == NULL)
            {
                break;
            }
        }
    }

    // If we have not fixed all vertices we check if some vertices share the
    // same coordinate
    if (not_fixed_verts.size() > 0)
    {
        // Collect verts with same coordinate.
        std::vector<std::pair<int, int> > same_coord(0);
        float tol = 1e-3;

        for (int i = 0; i < not_fixed_verts.size() - 1; i++)
        {
            for (int j = i + 1; j < not_fixed_verts.size(); j++)
            {
                if ((std::fabs(vertex[not_fixed_verts[i]].x -
                               vertex[not_fixed_verts[j]].x) < tol) &&
                    (std::fabs(vertex[not_fixed_verts[i]].y -
                               vertex[not_fixed_verts[j]].y) < tol) &&
                    (std::fabs(vertex[not_fixed_verts[i]].z -
                               vertex[not_fixed_verts[j]].z) < tol))
                {
                    same_coord.push_back(std::make_pair(not_fixed_verts[i],
                                                        not_fixed_verts[j]));
                    printf("Same coordinates: (%d, %d)\n", not_fixed_verts[i],
                           not_fixed_verts[j]);
                }
            }
        }
    }

    printf("New number of vertices: %d\n", new_number_of_vertices);

    // No vertices added
    if (new_number_of_vertices == num_vertices)
    {
        return;
    }

    // Add collected vertices and neighbors
    FLTVECT *new_vertices = (FLTVECT *)malloc(sizeof(FLTVECT) * new_number_of_vertices);

    // Counter for added new entities
    m = 0;

    for (n = 0; n < new_number_of_vertices; n++)
    {
        // If we are not adding any vertices
        if (n < num_vertices)
        {
            // Copy the old vertices
            new_vertices[n].x   = vertex[n].x;
            new_vertices[n].y   = vertex[n].y;
            new_vertices[n].z   = vertex[n].z;
            new_vertices[n].sel = vertex[n].sel;
            new_vertices[n].m   = vertex[n].m;
        }
        else
        {
            // Add new vertices
            new_vertices[n].x   = tmp_vertices[m].x;
            new_vertices[n].y   = tmp_vertices[m].y;
            new_vertices[n].z   = tmp_vertices[m].z;
            new_vertices[n].sel = tmp_vertices[m].sel;
            new_vertices[n].m   = tmp_vertices[m].m;

            m++;
        }
    }

    // Delete old data and reassign with new
    destroyNeighborlist();

    free(vertex);
    vertex       = new_vertices;
    num_vertices = new_number_of_vertices;
    closed       = true;

    // Debug
    // SurfaceMeshOld::writeOFF(surfmesh, "new_surfmesh.off");

    // Re-create neighbors
    createNeighborlist();
}

void SurfaceMeshOld::removeUnconnectedPatches(int minimal_number)
{
    // Re-create the neigborlist
    createNeighborlist();

    // Scratch data
    std::vector<int>  will_visit(0);
    std::vector<int>  patches_faces(num_faces, -1);
    std::vector<int>  patches_vertices(num_vertices, -1);
    std::vector<bool> visited_verts(num_vertices, false);
    std::vector<bool> visited_faces(num_faces, false);
    will_visit.reserve(num_vertices);

    int patch       = 0;
    int num_visited = 0;

    // Iterate over the different patches
    while (num_visited < num_vertices)
    {
        // Find next vertex which will be visited
        for (int i = 0; i < num_vertices; i++)
        {
            // Check if we have visited the vert
            if (visited_verts[i])
            {
                continue;
            }

            // Schedule the vert to be visited
            will_visit.push_back(i);
            break;
        }

        // Iterate over all neigboring vertices in a patch
        while (will_visit.size() > 0)
        {
            //  Bump number of visited verts
            num_visited += 1;

            // Vertex which we now visit
            const int visiting = will_visit.back();
            will_visit.pop_back();

            // printf("Patch: %d, Num visited: %d, Now visiting: %d, # to visit: %d\n",
            //       patch, num_visited, visiting, will_visit.size());

            // Store the patch
            patches_vertices[visiting] = patch;

            // Get all neigbors
            NPNT3 *tmp_ngr = neighbor_list[visiting];

            while (tmp_ngr != NULL)
            {
                if (!visited_verts[tmp_ngr->a])
                {
                    will_visit.push_back(tmp_ngr->a);
                    visited_verts[tmp_ngr->a] = true;
                }

                patches_faces[tmp_ngr->c] = patch;
                tmp_ngr                   = tmp_ngr->next;
            }
        }

        // Bump the patch number
        patch += 1;
    }

    printf("Num patches: %d\n", patch);

    // Collect the number faces in each patch
    std::vector<int> num_faces_per_patch(patch, 0);

    for (int i = 0; i < num_faces; i++)
    {
        num_faces_per_patch[patches_faces[i]] += 1;
    }

    for (int i = 0; i < patch; i++)
    {
        printf("Num faces in patch %d: %d\n", i, num_faces_per_patch[i]);
    }
}
