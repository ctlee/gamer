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

#include "biom.h"
#include "SurfaceMeshOld.h"
#include <limits>

/*
 * ***************************************************************************
 * Routine:  SurfaceMeshOld::correctNormals
 *
 * Author:   Zeyun Yu (zeyun.yu@gmail.com)
 *
 * Purpose:  Function to correct normals so all facet normals point outwards
 * ***************************************************************************
 */
void SurfaceMeshOld::correctNormals()
{
    int  n;
    int  i, j;
    int  a, b, c;
    int  a0, b0, c0;
    int *stack;
    unsigned char *visited;
    int   start, end;
    float x, y, z;
    float distance, max_dist;
    int   progress;
    float max_x, max_y, max_z;
    float ax, ay, az;
    float bx, by, bz;
    float cx, cy, cz;

    // get the max coordinates of the nodes
    max_x = -std::numeric_limits<float>::max();
    max_y = -std::numeric_limits<float>::max();
    max_z = -std::numeric_limits<float>::max();

    for (n = 0; n < num_vertices; n++)
    {
        x = vertex[n].x;
        y = vertex[n].y;
        z = vertex[n].z;

        if (x > max_x)
        {
            max_x = x;
        }

        if (y > max_y)
        {
            max_y = y;
        }

        if (z > max_z)
        {
            max_z = z;
        }
    }

    // printf("max coordinates: %f %f %f\n",max_x,max_y,max_z);
    max_x += 10.0f;
    max_y += 10.0f;
    max_z += 10.0f;

    // update the normal directions
    stack   = (int *)malloc(sizeof(int) * num_faces);
    visited = (unsigned char *)malloc(sizeof(unsigned char) * num_faces);

    for (i = 0; i < num_faces; i++)
    {
        visited[i] = 0;
    }
    start    = 0;
    progress = 0;

    while (1)
    {
        // find the unvisited face that is closest to (max_x,max_y,max_z)
        max_dist = 999999.0f;

        for (i = 0; i < num_faces; i++)
        {
            if (visited[i] == 0)
            {
                a        = face[i].a;
                b        = face[i].b;
                c        = face[i].c;
                x        = (vertex[a].x + vertex[b].x + vertex[c].x) / 3.0f;
                y        = (vertex[a].y + vertex[b].y + vertex[c].y) / 3.0f;
                z        = (vertex[a].z + vertex[b].z + vertex[c].z) / 3.0f;
                distance = sqrt((x - max_x) * (x - max_x) + (y - max_y) * (y - max_y) + (z - max_z) * (z - max_z));

                if (distance < max_dist)
                {
                    j        = i;
                    max_dist = distance;
                }
            }
        }

        if (max_dist == 999999.0f)
        {
            break;
        }

        a   = face[j].a;
        b   = face[j].b;
        c   = face[j].c;
        ax  = vertex[a].x;
        ay  = vertex[a].y;
        az  = vertex[a].z;
        bx  = vertex[b].x;
        by  = vertex[b].y;
        bz  = vertex[b].z;
        cx  = vertex[c].x;
        cy  = vertex[c].y;
        cz  = vertex[c].z;
        x   = (ax + bx + cx) / 3.0f;
        y   = (ay + by + cy) / 3.0f;
        z   = (az + bz + cz) / 3.0f;
        bx -= ax;
        by -= ay;
        bz -= az;
        cx -= ax;
        cy -= ay;
        cz -= az;
        ax  = by * cz - bz * cy;
        ay  = bz * cx - bx * cz;
        az  = bx * cy - by * cx;

        if (ax * (max_x - x) + ay * (max_y - y) + az * (max_z - z) < 0)
        {
            // switch the order of the last two vertices
            b         = face[j].b;
            face[j].b = face[j].c;
            face[j].c = b;
        }

        progress += start;
        stack[0]  = j;
        start     = 0;
        end       = 1;

        while (start < end)
        {
            if (((start) % 500) == 0)
            {
                printf("%2.2f%% done (%08d)\r", 100.0 * (float)(start + progress) / (float)(num_faces), start + progress);
                fflush(stdout);
            }

            i          = stack[start];
            visited[i] = 1;
            start++;
            a0 = face[i].a;
            b0 = face[i].b;
            c0 = face[i].c;

            for (n = 0; n < num_faces; n++)
            {
                if (visited[n] != 1)
                {
                    a = face[n].a;
                    b = face[n].b;
                    c = face[n].c;

                    if (visited[n] == 0)
                    {
                        if ((a == a0) && (b == b0))
                        {
                            j          = face[n].a;
                            face[n].a  = face[n].b;
                            face[n].b  = j;
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((b == a0) && (a == b0))
                        {
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((c == a0) && (a == b0))
                        {
                            j          = face[n].a;
                            face[n].a  = face[n].c;
                            face[n].c  = j;
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((a == a0) && (c == b0))
                        {
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((b == a0) && (c == b0))
                        {
                            j          = face[n].b;
                            face[n].b  = face[n].c;
                            face[n].c  = j;
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((c == a0) && (b == b0))
                        {
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }

                        else if ((b == a0) && (a == c0))
                        {
                            j          = face[n].a;
                            face[n].a  = face[n].b;
                            face[n].b  = j;
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((a == a0) && (b == c0))
                        {
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((a == a0) && (c == c0))
                        {
                            j          = face[n].a;
                            face[n].a  = face[n].c;
                            face[n].c  = j;
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((c == a0) && (a == c0))
                        {
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((c == a0) && (b == c0))
                        {
                            j          = face[n].b;
                            face[n].b  = face[n].c;
                            face[n].c  = j;
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((b == a0) && (c == c0))
                        {
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }

                        else if ((a == b0) && (b == c0))
                        {
                            j          = face[n].a;
                            face[n].a  = face[n].b;
                            face[n].b  = j;
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((b == b0) && (a == c0))
                        {
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((c == b0) && (a == c0))
                        {
                            j          = face[n].a;
                            face[n].a  = face[n].c;
                            face[n].c  = j;
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((a == b0) && (c == c0))
                        {
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((b == b0) && (c == c0))
                        {
                            j          = face[n].b;
                            face[n].b  = face[n].c;
                            face[n].c  = j;
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                        else if ((c == b0) && (b == c0))
                        {
                            visited[n] = 1;
                            stack[end] = n;
                            end++;
                        }
                    }
                }
            }
        }
    }
}
