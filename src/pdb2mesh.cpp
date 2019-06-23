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

#include "gamer/SurfaceMesh.h"
#include "gamer/PDBReader.h"
#include <vector>
#include <memory>
#include <cmath>

/// Namespace for all things gamer
namespace gamer
{


#define IndexVect1(i, j, k) ((k) * xdim1 * ydim1 + (j) * xdim1 + (i))
#define MaxVal            999999
#define MaxAtom           10

struct MOL_VERTEX {
    float         x;     // vertex coordinate
    float         y;
    float         z;
    unsigned char neigh; // first bit: +x; second bit: -x
                         // third bit: +y; fourth bit: -y
                         // fifth bit: +z; sixth  bit: -z
    unsigned short px;   // the corresponding index
    unsigned short py;
    unsigned short pz;
};

/** @brief Other data structure FLTVECT (float) */
struct FLTVECT {
    float x;   /**< @brief x-coordinate */
    float y;   /**< @brief y-coordinate */
    float z;   /**< @brief z-coordinate */
    int   m;   /**< @brief Marker */
    bool  sel; /**< @brief selection flag */
};


/** @brief Other data structure FLT2VECT (float) */
struct FLT2VECT {
    float x; /**< @brief x-coordinate */
    float y; /**< @brief y-coordinate */
};

/** @brief Other data structure INT4VECT (int) */
struct INT4VECT {
    int a; /**< @brief first integer */
    int b; /**< @brief second integer */
    int c; /**< @brief third integer */
    int d; /**< @brief fourth integer */
};

/** @brief Other data structure MinHeapS */
struct MinHeapS {
    unsigned short *x;    /**< @brief x-coordinate */
    unsigned short *y;    /**< @brief y-coordinate */
    unsigned short *z;    /**< @brief z-coordinate */
    int            *seed; /**< @brief seed */
    float          *dist; /**< @brief distance */
    int             size; /**< @brief size */
};

/** @brief Other data structure SEEDS */
struct SEEDS {
    float seedx;         /**< @brief x-coordinate */
    float seedy;         /**< @brief y-coordinate */
    float seedz;         /**< @brief z-coordinate */
    int   atom[MaxAtom]; /**< @brief atom array */
};

/** @brief Other data structure ATOM */
struct ATOM {
    float x;      /**< @brief x-coordinate */
    float y;      /**< @brief y-coordinate */
    float z;      /**< @brief z-coordinate */
    float radius; /**< @brief radius */
};



int GLOBAL_xdim, GLOBAL_ydim, GLOBAL_zdim;

// GRID variables
int *GLOBAL_segment_index;
int *GLOBAL_atom_index;

// Border variables
INT4VECT   *GLOBAL_quads;   int GLOBAL_vert_num;
MOL_VERTEX *GLOBAL_vertex;  int GLOBAL_quad_num;



#undef IndexVect
#define IndexVect(i, j, k) ( ( (k) * GLOBAL_ydim + (j) ) * GLOBAL_xdim + (i) )


int  CheckFaceCorner(float x,
                     float y,
                     float z);

void SetAtomIndex(int,
                  int,
                  int,
                  int);
int   ExtractSAS(int   atom_num,
                 ATOM *atom_list);
char  CheckManifold(int i,
                    int j,
                    int k);
float GetAngle(int a,
               int b,
               int c);


#define MaxDist    29999

int            xdim1, ydim1, zdim1;
MinHeapS      *min_heap;
SEEDS         *AllSeeds;
int           *heap_pointer;
ATOM          *atom_list;
float          threshold;
unsigned short min_x, min_y, min_z;
int            min_seed;
float          min_dist;


void GetMinimum(void);
void InsertHeap(int,
                int,
                int,
                float);
void UpdateHeap(int,
                int,
                int,
                float);
void    Marching(void);
FLTVECT FindSeed(float,
                 float,
                 float,
                 int);


void ExtractSES(MinHeapS *mheap, SEEDS *all_seeds, int *heappointer, int xd, int yd, int zd,
                int *atom_index, int atom_num, ATOM *atomlist, float thresh)
{
    int     i, j, k;
    int     m, n, l, num, c;
    int     index, index1;
    float   dist;
    FLTVECT seed;
    char    visited;


    xdim1 = xd;
    ydim1 = yd;
    zdim1 = zd;
    atom_list = atomlist;
    min_heap  = mheap;
    AllSeeds  = all_seeds;
    heap_pointer = heappointer;
    threshold = thresh;

    /* Initialize */
    index = 0;
    min_heap->size = 0;

    for (k = 0; k < zdim1; k++)
    {
        for (j = 0; j < ydim1; j++)
        {
            for (i = 0; i < xdim1; i++)
            {
                if (atom_index[IndexVect1(i, j, k)] < 0)
                {
                    for (num = 0; num < MaxAtom; num++)
                    {
                        AllSeeds[index].atom[num] = -1;
                    }
                    num = 0;

                    for (l = k - 1; l <= k + 1; l++)
                    {
                        for (n = j - 1; n <= j + 1; n++)
                        {
                            for (m = i - 1; m <= i + 1; m++)
                            {
                                if ((m == i) || (n == j) || (l == k))
                                {
                                    index1 = atom_index[IndexVect1(m, n, l)];

                                    if (index1 < 0)
                                    {
                                        index1  = -index1 - 1;
                                        visited = 0;

                                        for (c = 0; c < num; c++)
                                        {
                                            if (index1 == AllSeeds[index].atom[c])
                                            {
                                                visited = 1;
                                            }
                                        }

                                        if (visited == 0)
                                        {
                                            AllSeeds[index].atom[num] = index1;
                                            num++;

                                            if (num == MaxAtom)
                                            {
                                                num--;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    seed = FindSeed(i, j, k, index);
                    AllSeeds[index].seedx = seed.x;
                    AllSeeds[index].seedy = seed.y;
                    AllSeeds[index].seedz = seed.z;
                    dist = (seed.x - i) * (seed.x - i) + (seed.y - j) * (seed.y - j) + (seed.z - k) * (seed.z - k);
                    min_seed = index;
                    InsertHeap(i, j, k, dist);

                    index++;
                }
                else if (atom_index[IndexVect1(i, j, k)] > 0)
                {
                    heap_pointer[IndexVect1(i, j, k)] = MaxVal;
                }
                else
                {
                    heap_pointer[IndexVect1(i, j, k)] = -11;
                }
            }
        }
    }


    /* Fast Marching Method */
    while (1)
    {
        GetMinimum();

        if (min_dist >= MaxDist - 0.001)
        {
            break;
        }

        Marching();
    }
}

void GetMinimum(void)
{
    int   pointer, left, right;
    float dist;

    min_x = min_heap->x[0];
    min_y = min_heap->y[0];
    min_z = min_heap->z[0];
    min_seed = min_heap->seed[0];
    min_dist = min_heap->dist[0];


    if (min_dist == MaxDist)
    {
        return;
    }

    heap_pointer[IndexVect1(min_heap->x[0], min_heap->y[0], min_heap->z[0])] = -3;

    min_heap->size--;
    dist = min_heap->dist[min_heap->size];

    pointer = 1;

    while (pointer <= min_heap->size / 2)
    {
        left  = 2 * pointer;
        right = 2 * pointer + 1;

        if ((min_heap->dist[left - 1] <= min_heap->dist[right - 1]) && (min_heap->dist[left - 1] < dist))
        {
            min_heap->x[pointer - 1] = min_heap->x[left - 1];
            min_heap->y[pointer - 1] = min_heap->y[left - 1];
            min_heap->z[pointer - 1] = min_heap->z[left - 1];
            min_heap->seed[pointer - 1] = min_heap->seed[left - 1];
            min_heap->dist[pointer - 1] = min_heap->dist[left - 1];
            heap_pointer[IndexVect1(min_heap->x[pointer - 1], min_heap->y[pointer - 1], min_heap->z[pointer - 1])] = pointer - 1;
            pointer = left;
        }
        else if ((min_heap->dist[left - 1] > min_heap->dist[right - 1]) && (min_heap->dist[right - 1] < dist))
        {
            min_heap->x[pointer - 1] = min_heap->x[right - 1];
            min_heap->y[pointer - 1] = min_heap->y[right - 1];
            min_heap->z[pointer - 1] = min_heap->z[right - 1];
            min_heap->seed[pointer - 1] = min_heap->seed[right - 1];
            min_heap->dist[pointer - 1] = min_heap->dist[right - 1];
            heap_pointer[IndexVect1(min_heap->x[pointer - 1], min_heap->y[pointer - 1], min_heap->z[pointer - 1])] = pointer - 1;
            pointer = right;
        }
        else
        {
            break;
        }
    }

    min_heap->x[pointer - 1] = min_heap->x[min_heap->size];
    min_heap->y[pointer - 1] = min_heap->y[min_heap->size];
    min_heap->z[pointer - 1] = min_heap->z[min_heap->size];
    min_heap->seed[pointer - 1] = min_heap->seed[min_heap->size];
    min_heap->dist[pointer - 1] = dist;
    heap_pointer[IndexVect1(min_heap->x[min_heap->size], min_heap->y[min_heap->size], min_heap->z[min_heap->size])] = pointer - 1;
}

void InsertHeap(int x, int y, int z, float dist)
{
    int parent = 0;


    min_heap->size++;
    int pointer = min_heap->size;

    while (pointer > 1)
    {
        if (pointer % 2 == 0)
        {
            parent = pointer / 2;
        }
        else if (pointer % 2 == 1)
        {
            parent = (pointer - 1) / 2;
        }

        if (dist < min_heap->dist[parent - 1])
        {
            min_heap->x[pointer - 1] = min_heap->x[parent - 1];
            min_heap->y[pointer - 1] = min_heap->y[parent - 1];
            min_heap->z[pointer - 1] = min_heap->z[parent - 1];
            min_heap->seed[pointer - 1] = min_heap->seed[parent - 1];
            min_heap->dist[pointer - 1] = min_heap->dist[parent - 1];

            int index = IndexVect1(min_heap->x[pointer - 1], min_heap->y[pointer - 1], min_heap->z[pointer - 1]);
            heap_pointer[index] = pointer - 1;

            pointer = parent;
        }
        else
        {
            break;
        }
    }
    min_heap->x[pointer - 1] = x;
    min_heap->y[pointer - 1] = y;
    min_heap->z[pointer - 1] = z;
    min_heap->seed[pointer - 1] = min_seed;
    min_heap->dist[pointer - 1] = dist;

    heap_pointer[IndexVect1(x, y, z)] = pointer - 1;
}

void UpdateHeap(int x, int y, int z, float dist)
{
    int  parent = 0;
    int  left, right;
    char up = 0;


    int pointer = heap_pointer[IndexVect1(x, y, z)] + 1;

    // checking the upper elements
    while (pointer > 1)
    {
        if (pointer % 2 == 0)
        {
            parent = pointer / 2;
        }
        else if (pointer % 2 == 1)
        {
            parent = (pointer - 1) / 2;
        }

        if (dist < min_heap->dist[parent - 1])
        {
            min_heap->x[pointer - 1] = min_heap->x[parent - 1];
            min_heap->y[pointer - 1] = min_heap->y[parent - 1];
            min_heap->z[pointer - 1] = min_heap->z[parent - 1];
            min_heap->seed[pointer - 1] = min_heap->seed[parent - 1];
            min_heap->dist[pointer - 1] = min_heap->dist[parent - 1];

            int index = IndexVect1(min_heap->x[pointer - 1], min_heap->y[pointer - 1], min_heap->z[pointer - 1]);
            heap_pointer[index] = pointer - 1;

            pointer = parent;
        }
        else
        {
            break;
        }
    }

    if (up == 0)
    {
        // checking the lower elements
        while (pointer <= min_heap->size / 2)
        {
            left  = 2 * pointer;
            right = 2 * pointer + 1;

            if ((min_heap->dist[left - 1] <= min_heap->dist[right - 1]) && (min_heap->dist[left - 1] < dist))
            {
                min_heap->x[pointer - 1] = min_heap->x[left - 1];
                min_heap->y[pointer - 1] = min_heap->y[left - 1];
                min_heap->z[pointer - 1] = min_heap->z[left - 1];
                min_heap->seed[pointer - 1] = min_heap->seed[left - 1];
                min_heap->dist[pointer - 1] = min_heap->dist[left - 1];
                heap_pointer[IndexVect1(min_heap->x[pointer - 1], min_heap->y[pointer - 1], min_heap->z[pointer - 1])] = pointer - 1;
                pointer = left;
            }
            else if ((min_heap->dist[left - 1] > min_heap->dist[right - 1]) && (min_heap->dist[right - 1] < dist))
            {
                min_heap->x[pointer - 1] = min_heap->x[right - 1];
                min_heap->y[pointer - 1] = min_heap->y[right - 1];
                min_heap->z[pointer - 1] = min_heap->z[right - 1];
                min_heap->seed[pointer - 1] = min_heap->seed[right - 1];
                min_heap->dist[pointer - 1] = min_heap->dist[right - 1];
                heap_pointer[IndexVect1(min_heap->x[pointer - 1], min_heap->y[pointer - 1], min_heap->z[pointer - 1])] = pointer - 1;
                pointer = right;
            }
            else
            {
                break;
            }
        }
    }


    min_heap->x[pointer - 1] = x;
    min_heap->y[pointer - 1] = y;
    min_heap->z[pointer - 1] = z;
    min_heap->seed[pointer - 1] = min_seed;
    min_heap->dist[pointer - 1] = dist;

    heap_pointer[IndexVect1(x, y, z)] = pointer - 1;
}

void Marching(void)
{
    int   tempt_x, tempt_y, tempt_z;
    float dt, dist;
    int   neighbor, seed;
    char  boundary;
    float min_seedx, min_seedy, min_seedz;
    float seedx, seedy, seedz;


    min_seedx = AllSeeds[min_seed].seedx;
    min_seedy = AllSeeds[min_seed].seedy;
    min_seedz = AllSeeds[min_seed].seedz;

    boundary = 0;


    /* ============================ */
    tempt_x = std::max(min_x - 1, 0);
    tempt_y = min_y;
    tempt_z = min_z;

    if (heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)] == MaxVal)
    {
        dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) * (tempt_z - min_seedz);

        if (dist <= threshold)
        {
            InsertHeap(tempt_x, tempt_y, tempt_z, dist);
        }
        else
        {
            boundary = 1;
        }
    }
    else if (heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)] > -1)
    {
        neighbor = heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)];
        dt = min_heap->dist[neighbor];

        if (dt < MaxDist)
        {
            dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) *
                (tempt_z - min_seedz);

            if (dist < dt)
            {
                UpdateHeap(tempt_x, tempt_y, tempt_z, dist);
            }
        }
        else
        {
            dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) *
                (tempt_z - min_seedz);
            seed  = min_heap->seed[neighbor];
            seedx = AllSeeds[seed].seedx;
            seedy = AllSeeds[seed].seedy;
            seedz = AllSeeds[seed].seedz;

            if (dist < (tempt_x - seedx) * (tempt_x - seedx) + (tempt_y - seedy) * (tempt_y - seedy) + (tempt_z - seedz) * (tempt_z - seedz))
            {
                UpdateHeap(tempt_x, tempt_y, tempt_z, MaxDist);
            }
        }
    }

    // else if (heap_pointer[IndexVect1(tempt_x,tempt_y,tempt_z)] > -11) {


    /* ============================ */
    tempt_x = std::min(min_x + 1, xdim1 - 1);
    tempt_y = min_y;
    tempt_z = min_z;

    if (heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)] == MaxVal)
    {
        dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) * (tempt_z - min_seedz);

        if (dist <= threshold)
        {
            InsertHeap(tempt_x, tempt_y, tempt_z, dist);
        }
        else
        {
            boundary = 1;
        }
    }
    else if (heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)] > -1)
    {
        neighbor = heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)];
        dt = min_heap->dist[neighbor];

        if (dt < MaxDist)
        {
            dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) *
                (tempt_z - min_seedz);

            if (dist < dt)
            {
                UpdateHeap(tempt_x, tempt_y, tempt_z, dist);
            }
        }
        else
        {
            dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) *
                (tempt_z - min_seedz);
            seed  = min_heap->seed[neighbor];
            seedx = AllSeeds[seed].seedx;
            seedy = AllSeeds[seed].seedy;
            seedz = AllSeeds[seed].seedz;

            if (dist < (tempt_x - seedx) * (tempt_x - seedx) + (tempt_y - seedy) * (tempt_y - seedy) + (tempt_z - seedz) * (tempt_z - seedz))
            {
                UpdateHeap(tempt_x, tempt_y, tempt_z, MaxDist);
            }
        }
    }


    /* ============================ */
    tempt_y = std::max(min_y - 1, 0);
    tempt_x = min_x;
    tempt_z = min_z;

    if (heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)] == MaxVal)
    {
        dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) * (tempt_z - min_seedz);

        if (dist <= threshold)
        {
            InsertHeap(tempt_x, tempt_y, tempt_z, dist);
        }
        else
        {
            boundary = 1;
        }
    }
    else if (heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)] > -1)
    {
        neighbor = heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)];
        dt = min_heap->dist[neighbor];

        if (dt < MaxDist)
        {
            dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) *
                (tempt_z - min_seedz);

            if (dist < dt)
            {
                UpdateHeap(tempt_x, tempt_y, tempt_z, dist);
            }
        }
        else
        {
            dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) *
                (tempt_z - min_seedz);
            seed  = min_heap->seed[neighbor];
            seedx = AllSeeds[seed].seedx;
            seedy = AllSeeds[seed].seedy;
            seedz = AllSeeds[seed].seedz;

            if (dist < (tempt_x - seedx) * (tempt_x - seedx) + (tempt_y - seedy) * (tempt_y - seedy) + (tempt_z - seedz) * (tempt_z - seedz))
            {
                UpdateHeap(tempt_x, tempt_y, tempt_z, MaxDist);
            }
        }
    }


    /* ============================ */
    tempt_y = std::min(min_y + 1, ydim1 - 1);
    tempt_x = min_x;
    tempt_z = min_z;

    if (heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)] == MaxVal)
    {
        dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) * (tempt_z - min_seedz);

        if (dist <= threshold)
        {
            InsertHeap(tempt_x, tempt_y, tempt_z, dist);
        }
        else
        {
            boundary = 1;
        }
    }
    else if (heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)] > -1)
    {
        neighbor = heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)];
        dt = min_heap->dist[neighbor];

        if (dt < MaxDist)
        {
            dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) *
                (tempt_z - min_seedz);

            if (dist < dt)
            {
                UpdateHeap(tempt_x, tempt_y, tempt_z, dist);
            }
        }
        else
        {
            dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) *
                (tempt_z - min_seedz);
            seed  = min_heap->seed[neighbor];
            seedx = AllSeeds[seed].seedx;
            seedy = AllSeeds[seed].seedy;
            seedz = AllSeeds[seed].seedz;

            if (dist < (tempt_x - seedx) * (tempt_x - seedx) + (tempt_y - seedy) * (tempt_y - seedy) + (tempt_z - seedz) * (tempt_z - seedz))
            {
                UpdateHeap(tempt_x, tempt_y, tempt_z, MaxDist);
            }
        }
    }


    /* ============================ */
    tempt_z = std::max(min_z - 1, 0);
    tempt_x = min_x;
    tempt_y = min_y;

    if (heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)] == MaxVal)
    {
        dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) * (tempt_z - min_seedz);

        if (dist <= threshold)
        {
            InsertHeap(tempt_x, tempt_y, tempt_z, dist);
        }
        else
        {
            boundary = 1;
        }
    }
    else if (heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)] > -1)
    {
        neighbor = heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)];
        dt = min_heap->dist[neighbor];

        if (dt < MaxDist)
        {
            dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) *
                (tempt_z - min_seedz);

            if (dist < dt)
            {
                UpdateHeap(tempt_x, tempt_y, tempt_z, dist);
            }
        }
        else
        {
            dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) *
                (tempt_z - min_seedz);
            seed  = min_heap->seed[neighbor];
            seedx = AllSeeds[seed].seedx;
            seedy = AllSeeds[seed].seedy;
            seedz = AllSeeds[seed].seedz;

            if (dist < (tempt_x - seedx) * (tempt_x - seedx) + (tempt_y - seedy) * (tempt_y - seedy) + (tempt_z - seedz) * (tempt_z - seedz))
            {
                UpdateHeap(tempt_x, tempt_y, tempt_z, MaxDist);
            }
        }
    }


    /* ============================ */
    tempt_z = std::min(min_z + 1, zdim1 - 1);
    tempt_x = min_x;
    tempt_y = min_y;

    if (heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)] == MaxVal)
    {
        dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) * (tempt_z - min_seedz);

        if (dist <= threshold)
        {
            InsertHeap(tempt_x, tempt_y, tempt_z, dist);
        }
        else
        {
            boundary = 1;
        }
    }
    else if (heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)] > -1)
    {
        neighbor = heap_pointer[IndexVect1(tempt_x, tempt_y, tempt_z)];
        dt = min_heap->dist[neighbor];

        if (dt < MaxDist)
        {
            dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) *
                (tempt_z - min_seedz);

            if (dist < dt)
            {
                UpdateHeap(tempt_x, tempt_y, tempt_z, dist);
            }
        }
        else
        {
            dist = (tempt_x - min_seedx) * (tempt_x - min_seedx) + (tempt_y - min_seedy) * (tempt_y - min_seedy) + (tempt_z - min_seedz) *
                (tempt_z - min_seedz);
            seed  = min_heap->seed[neighbor];
            seedx = AllSeeds[seed].seedx;
            seedy = AllSeeds[seed].seedy;
            seedz = AllSeeds[seed].seedz;

            if (dist < (tempt_x - seedx) * (tempt_x - seedx) + (tempt_y - seedy) * (tempt_y - seedy) + (tempt_z - seedz) * (tempt_z - seedz))
            {
                UpdateHeap(tempt_x, tempt_y, tempt_z, MaxDist);
            }
        }
    }


    if (boundary)
    {
        InsertHeap(min_x, min_y, min_z, MaxDist);
    }
}

FLTVECT FindSeed(float x, float y, float z, int index)
{
    double  cx1, cy1, cz1;
    double  cx2, cy2, cz2;
    double  cx3, cy3, cz3;
    double  dist, radius1, radius2, radius3;
    double  cos_alpha;
    double  ax, ay, az;
    double  bx, by, bz;
    double  dx, dy, dz;
    double  hx, hy, hz;
    int     atom1, atom2, atom3;
    int     num, total;
    FLTVECT tmp;


    // AllSeeds[index].atom[0] is always >= 0
    atom1 = AllSeeds[index].atom[0];
    atom2 = AllSeeds[index].atom[1];
    atom3 = AllSeeds[index].atom[2];

    // contain only one atom
    if (atom2 < 0)
    {
        radius1 = atom_list[atom1].radius;
        cx1  = atom_list[atom1].x;
        cy1  = atom_list[atom1].y;
        cz1  = atom_list[atom1].z;
        dx   = x - cx1;
        dy   = y - cy1;
        dz   = z - cz1;
        dist = sqrt(dx * dx + dy * dy + dz * dz);

        tmp.x = cx1 + radius1 * dx / dist;
        tmp.y = cy1 + radius1 * dy / dist;
        tmp.z = cz1 + radius1 * dz / dist;
    }

    // contain only two atoms
    else if (atom3 < 0)
    {
        radius1 = atom_list[atom1].radius;
        cx1 = atom_list[atom1].x;
        cy1 = atom_list[atom1].y;
        cz1 = atom_list[atom1].z;
        radius2 = atom_list[atom2].radius;
        cx2  = atom_list[atom2].x;
        cy2  = atom_list[atom2].y;
        cz2  = atom_list[atom2].z;
        dist = sqrt((cx2 - cx1) * (cx2 - cx1) + (cy2 - cy1) * (cy2 - cy1) + (cz2 - cz1) * (cz2 - cz1));
        cos_alpha = (radius1 * radius1 + dist * dist - radius2 * radius2) / (2.0 * radius1 * dist);

        ax   = (cx2 - cx1) / dist;
        ay   = (cy2 - cy1) / dist;
        az   = (cz2 - cz1) / dist;
        bx   = x - cx1;
        by   = y - cy1;
        bz   = z - cz1;
        dx   = ay * bz - az * by;
        dy   = az * bx - ax * bz;
        dz   = ax * by - ay * bx;
        hx   = dy * az - dz * ay;
        hy   = dz * ax - dx * az;
        hz   = dx * ay - dy * ax;
        dist = sqrt(hx * hx + hy * hy + hz * hz);
        hx  /= dist;
        hy  /= dist;
        hz  /= dist;
        dist = radius1 * sqrt(1.0 - cos_alpha * cos_alpha);

        radius1 *= cos_alpha;
        tmp.x    = radius1 * ax + cx1 + dist * hx;
        tmp.y    = radius1 * ay + cy1 + dist * hy;
        tmp.z    = radius1 * az + cz1 + dist * hz;
    }

    // contain three or more atoms
    else
    {
        hx = 0;
        hy = 0;
        hz = 0;
        total = 0;
        num   = 2;

        while (num < MaxAtom)
        {
            if (AllSeeds[index].atom[num] < 0)
            {
                break;
            }
            atom1 = AllSeeds[index].atom[num - 2];
            atom2 = AllSeeds[index].atom[num - 1];
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

            while (1)
            {
                // projection to the first atom sphere
                dx   = ax - cx1;
                dy   = ay - cy1;
                dz   = az - cz1;
                dist = radius1 / sqrt(dx * dx + dy * dy + dz * dz);
                ax   = cx1 + dist * dx;
                ay   = cy1 + dist * dy;
                az   = cz1 + dist * dz;

                // projection to the second atom sphere
                dx   = ax - cx2;
                dy   = ay - cy2;
                dz   = az - cz2;
                dist = radius2 / sqrt(dx * dx + dy * dy + dz * dz);
                ax   = cx2 + dist * dx;
                ay   = cy2 + dist * dy;
                az   = cz2 + dist * dz;

                // projection to the third atom sphere
                dx   = ax - cx3;
                dy   = ay - cy3;
                dz   = az - cz3;
                dist = radius3 / sqrt(dx * dx + dy * dy + dz * dz);
                ax   = cx3 + dist * dx;
                ay   = cy3 + dist * dy;
                az   = cz3 + dist * dz;

                dist = sqrt((ax - bx) * (ax - bx) + (ay - by) * (ay - by) + (az - bz) * (az - bz));

                if (dist < 0.001)
                {
                    hx += ax;
                    hy += ay;
                    hz += az;
                    total++;
                    break;
                }
                else
                {
                    bx = ax;
                    by = ay;
                    bz = az;
                }
            }
            num++;
        }
        tmp.x = hx / (float)total;
        tmp.y = hy / (float)total;
        tmp.z = hz / (float)total;
    }

    return tmp;
}

/*
 * ***************************************************************************
 * Routine:  getMinMax    < ... >
 *
 * Authors:   Zeyun Yu (zeyun.yu@gmail.com)
 *            John Moody (brogan@gmail.com)
 *
 * Purpose:  Calculate the minimum and maximum range of the blurred atoms
 * ***************************************************************************
 */
template <typename Iterator>
void getMinMax(Iterator begin,
               Iterator end,
               float    min[3],
               float    max[3])
{
    float maxRad = 0.0;
    float tempRad;

    min[0] = min[1] = min[2] = std::numeric_limits<float>::infinity();
    max[0] = max[1] = max[2] = -std::numeric_limits<float>::infinity();
    maxRad = 0;

    for (auto curr = begin; curr != end; ++curr)
    {
        float x = curr->x;
        float y = curr->y;
        float z = curr->z;

        min[0] = std::min(min[0], x);
        max[0] = std::max(max[0], x);

        min[1] = std::min(min[1], y);
        max[1] = std::max(max[1], y);

        min[2] = std::min(min[2], z);
        max[2] = std::max(max[2], z);

        // Definition of surface boundary
        tempRad = curr->radius * sqrt(1.0 + log(pdbreader_detail::EPSILON) / BLOBBYNESS);

        maxRad = std::max(maxRad, tempRad);
    }

    for (int i = 0; i < 3; i++)
    {
        min[i] -= maxRad;
        max[i] += maxRad;
    }
}

std::unique_ptr<SurfaceMesh> readPDB_molsurf(const std::string &input_name)
{
    int               i, j, k;
    int               a, b, c, d;
    float             orig[3], span[3];
    int               dim[3];
    int               m, n, l, num;
    double            threshold;
    double            nx, ny, nz;
    int               xydim, xyzdim;
    clock_t           begin, finish;
    std::vector<Atom> atoms;
    std::vector<ATOM> atom_list;
    float             min[3], max[3];
    SEEDS            *AllSeeds; // Border variable
    MinHeapS         *min_heap;

    // Read in the PDB file
    readPDB(input_name, std::back_inserter(atoms));

    for (auto atom : atoms)
    {

        ATOM new_atom;
        new_atom.x = atom.pos[0];
        new_atom.y = atom.pos[1];
        new_atom.z = atom.pos[2];
        new_atom.radius = atom.radius;

        atom_list.push_back(new_atom);
    }

    getMinMax(atom_list.begin(), atom_list.end(), min, max);

    GLOBAL_xdim = (int)(((max[0] - min[0]) + 1) * DIM_SCALE);
    GLOBAL_ydim = (int)(((max[1] - min[1]) + 1) * DIM_SCALE);
    GLOBAL_zdim = (int)(((max[2] - min[2]) + 1) * DIM_SCALE);
    xydim  = GLOBAL_xdim * GLOBAL_ydim;
    xyzdim = xydim * GLOBAL_zdim;

    // printf("dimension: %d X %d X %d\n",xdim,ydim,zdim);

    GLOBAL_atom_index = (int *)malloc(sizeof(int) * xyzdim);
    GLOBAL_segment_index = (int *)malloc(sizeof(int) * xyzdim);

    for (k = 0; k < xyzdim; k++)
    {
        GLOBAL_atom_index[k] = 0;
    }

    orig[0] = min[0];
    orig[1] = min[1];
    orig[2] = min[2];
    span[0] = (max[0] - min[0]) / (double)(GLOBAL_xdim - 1);
    span[1] = (max[1] - min[1]) / (double)(GLOBAL_ydim - 1);
    span[2] = (max[2] - min[2]) / (double)(GLOBAL_zdim - 1);
    dim[0]  = GLOBAL_xdim;
    dim[1]  = GLOBAL_ydim;
    dim[2]  = GLOBAL_zdim;

    for (m = 0; m < atom_list.size(); m++)
    {
        atom_list[m].x = (atom_list[m].x - orig[0]) / span[0];
        atom_list[m].y = (atom_list[m].y - orig[1]) / span[1];
        atom_list[m].z = (atom_list[m].z - orig[2]) / span[2];
        atom_list[m].radius = (atom_list[m].radius + 1.5) / ((span[0] + span[1] + span[2]) / 3.0);
    }

    begin  = clock();
    num    = ExtractSAS(atom_list.size(), atom_list.data());
    finish = clock();

    // for(int q=0; q < GLOBAL_xdim; ++q){
    //     for(int r=0; r < GLOBAL_ydim; ++r){
    //         for(int s=0; s < GLOBAL_zdim; ++s){
    //             std::cout << GLOBAL_atom_index[IndexVect(q,r,s)] <<
    // std::endl;
    //         }
    //     }
    // }

    // printf("   Extract SAS voxels: CPU Time = %f seconds
    // \n",(double)(finish-begin)/CLOCKS_PER_SEC);
    // printf("   Number of boundary voxels: %d\n\n",num);


    begin = clock();
    threshold   = 1.5 / ((span[0] + span[1] + span[2]) / 3.0);
    min_heap    = (MinHeapS *)malloc(sizeof(MinHeapS));
    min_heap->x = (unsigned short *)malloc(sizeof(unsigned short) * num * 3);
    min_heap->y = (unsigned short *)malloc(sizeof(unsigned short) * num * 3);
    min_heap->z = (unsigned short *)malloc(sizeof(unsigned short) * num * 3);
    min_heap->seed = (int *)malloc(sizeof(int) * num * 3);
    min_heap->dist = (float *)malloc(sizeof(float) * num * 3);
    AllSeeds = (SEEDS *)malloc(sizeof(SEEDS) * num);
    ExtractSES(min_heap,
               AllSeeds,
               GLOBAL_segment_index,
               GLOBAL_xdim,
               GLOBAL_ydim,
               GLOBAL_zdim,
               GLOBAL_atom_index,
               atom_list.size(),
               atom_list.data(),
               threshold * threshold);
    finish = clock();

    // printf("   Extract SES voxels: CPU Time = %f seconds
    // \n\n",(double)(finish-begin)/CLOCKS_PER_SEC);


    // detect and fix non-manifolds !
    while (1)
    {
        b = 0;

        int index = 0;

        for (k = 0; k < GLOBAL_zdim; k++)
        {
            for (j = 0; j < GLOBAL_ydim; j++)
            {
                for (i = 0; i < GLOBAL_xdim; i++, ++index)
                {
                    if (GLOBAL_segment_index[index] == MaxVal)
                    {
                        if (!CheckManifold(i, j, k)) // non-manifold occurs
                        {
                            GLOBAL_segment_index[index] = 0;
                            min_heap->x[min_heap->size] = i;
                            min_heap->y[min_heap->size] = j;
                            min_heap->z[min_heap->size] = k;
                            min_heap->size++;
                            b++;
                        }
                    }
                }
            }
        }

        // printf("non-manifold = %d \n",b);
        if (b == 0)
        {
            break;
        }
    }


    // generate the surface mesh
    GLOBAL_vertex = (MOL_VERTEX *)malloc(sizeof(MOL_VERTEX) * num * 8);
    GLOBAL_quads  = (INT4VECT *)malloc(sizeof(INT4VECT) * num * 6);

    for (k = 0; k < num * 8; k++)
    {
        GLOBAL_vertex[k].neigh = 0;
    }

    for (k = 0; k < xyzdim; k++)
    {
        GLOBAL_atom_index[k] = -1;
    }
    begin = clock();
    GLOBAL_vert_num = 0;
    GLOBAL_quad_num = 0;

    for (num = 0; num < min_heap->size; num++)
    {
        i = min_heap->x[num];
        j = min_heap->y[num];
        k = min_heap->z[num];

        // back face
        if (GLOBAL_segment_index[IndexVect(i - 1, j, k)] == MaxVal)
        {
            a = CheckFaceCorner(i - 0.5, j - 0.5, k - 0.5);
            GLOBAL_vertex[a].neigh |= 40; // +y and +z
            b = CheckFaceCorner(i - 0.5, j - 0.5, k + 0.5);
            GLOBAL_vertex[b].neigh |= 24; // +y and -z
            c = CheckFaceCorner(i - 0.5, j + 0.5, k + 0.5);
            GLOBAL_vertex[c].neigh |= 20; // -y and -z
            d = CheckFaceCorner(i - 0.5, j + 0.5, k - 0.5);
            GLOBAL_vertex[d].neigh |= 36; // -y and +z

            GLOBAL_quads[GLOBAL_quad_num].a = a;
            GLOBAL_quads[GLOBAL_quad_num].b = b;
            GLOBAL_quads[GLOBAL_quad_num].c = c;
            GLOBAL_quads[GLOBAL_quad_num].d = d;
            GLOBAL_quad_num++;
        }

        // front face
        if (GLOBAL_segment_index[IndexVect(i + 1, j, k)] == MaxVal)
        {
            a = CheckFaceCorner(i + 0.5, j - 0.5, k - 0.5);
            GLOBAL_vertex[a].neigh |= 40; // +y and +z
            b = CheckFaceCorner(i + 0.5, j + 0.5, k - 0.5);
            GLOBAL_vertex[b].neigh |= 36; // -y and +z
            c = CheckFaceCorner(i + 0.5, j + 0.5, k + 0.5);
            GLOBAL_vertex[c].neigh |= 20; // -y and -z
            d = CheckFaceCorner(i + 0.5, j - 0.5, k + 0.5);
            GLOBAL_vertex[d].neigh |= 24; // +y and -z

            GLOBAL_quads[GLOBAL_quad_num].a = a;
            GLOBAL_quads[GLOBAL_quad_num].b = b;
            GLOBAL_quads[GLOBAL_quad_num].c = c;
            GLOBAL_quads[GLOBAL_quad_num].d = d;
            GLOBAL_quad_num++;
        }

        // left face
        if (GLOBAL_segment_index[IndexVect(i, j - 1, k)] == MaxVal)
        {
            a = CheckFaceCorner(i + 0.5, j - 0.5, k - 0.5);
            GLOBAL_vertex[a].neigh |= 33; // -x and +z
            b = CheckFaceCorner(i + 0.5, j - 0.5, k + 0.5);
            GLOBAL_vertex[b].neigh |= 17; // -x and -z
            c = CheckFaceCorner(i - 0.5, j - 0.5, k + 0.5);
            GLOBAL_vertex[c].neigh |= 18; // +x and -z
            d = CheckFaceCorner(i - 0.5, j - 0.5, k - 0.5);
            GLOBAL_vertex[d].neigh |= 34; // +x and +z

            GLOBAL_quads[GLOBAL_quad_num].a = a;
            GLOBAL_quads[GLOBAL_quad_num].b = b;
            GLOBAL_quads[GLOBAL_quad_num].c = c;
            GLOBAL_quads[GLOBAL_quad_num].d = d;
            GLOBAL_quad_num++;
        }

        // right face
        if (GLOBAL_segment_index[IndexVect(i, j + 1, k)] == MaxVal)
        {
            a = CheckFaceCorner(i + 0.5, j + 0.5, k - 0.5);
            GLOBAL_vertex[a].neigh |= 33; // -x and +z
            b = CheckFaceCorner(i - 0.5, j + 0.5, k - 0.5);
            GLOBAL_vertex[b].neigh |= 34; // +x and +z
            c = CheckFaceCorner(i - 0.5, j + 0.5, k + 0.5);
            GLOBAL_vertex[c].neigh |= 18; // +x and -z
            d = CheckFaceCorner(i + 0.5, j + 0.5, k + 0.5);
            GLOBAL_vertex[d].neigh |= 17; // -x and -z

            GLOBAL_quads[GLOBAL_quad_num].a = a;
            GLOBAL_quads[GLOBAL_quad_num].b = b;
            GLOBAL_quads[GLOBAL_quad_num].c = c;
            GLOBAL_quads[GLOBAL_quad_num].d = d;
            GLOBAL_quad_num++;
        }

        // bottom face
        if (GLOBAL_segment_index[IndexVect(i, j, k - 1)] == MaxVal)
        {
            a = CheckFaceCorner(i + 0.5, j - 0.5, k - 0.5);
            GLOBAL_vertex[a].neigh |= 9;  // -x and +y
            b = CheckFaceCorner(i - 0.5, j - 0.5, k - 0.5);
            GLOBAL_vertex[b].neigh |= 10; // +x and +y
            c = CheckFaceCorner(i - 0.5, j + 0.5, k - 0.5);
            GLOBAL_vertex[c].neigh |= 6;  // +x and -y
            d = CheckFaceCorner(i + 0.5, j + 0.5, k - 0.5);
            GLOBAL_vertex[d].neigh |= 5;  // -x and -y

            GLOBAL_quads[GLOBAL_quad_num].a = a;
            GLOBAL_quads[GLOBAL_quad_num].b = b;
            GLOBAL_quads[GLOBAL_quad_num].c = c;
            GLOBAL_quads[GLOBAL_quad_num].d = d;
            GLOBAL_quad_num++;
        }

        // top face
        if (GLOBAL_segment_index[IndexVect(i, j, k + 1)] == MaxVal)
        {
            a = CheckFaceCorner(i + 0.5, j - 0.5, k + 0.5);
            GLOBAL_vertex[a].neigh |= 9;  // -x and +y
            b = CheckFaceCorner(i + 0.5, j + 0.5, k + 0.5);
            GLOBAL_vertex[b].neigh |= 5;  // -x and -y
            c = CheckFaceCorner(i - 0.5, j + 0.5, k + 0.5);
            GLOBAL_vertex[c].neigh |= 6;  // +x and -y
            d = CheckFaceCorner(i - 0.5, j - 0.5, k + 0.5);
            GLOBAL_vertex[d].neigh |= 10; // +x and +y

            GLOBAL_quads[GLOBAL_quad_num].a = a;
            GLOBAL_quads[GLOBAL_quad_num].b = b;
            GLOBAL_quads[GLOBAL_quad_num].c = c;
            GLOBAL_quads[GLOBAL_quad_num].d = d;
            GLOBAL_quad_num++;
        }
    }
    finish = clock();

    // printf("   Generate quad meshes: CPU Time = %f seconds
    // \n",(double)(finish-begin)/CLOCKS_PER_SEC);
    // printf("   vert-num : %d -- quad-num: %d \n\n",vert_num,GLOBAL_quad_num);


    // Smooth the mesh
    begin = clock();
    unsigned char neighbor;

    for (num = 0; num < 3; num++)
    {
        for (n = 0; n < GLOBAL_vert_num; n++)
        {
            nx = 0;
            ny = 0;
            nz = 0;
            m  = 0;
            neighbor = GLOBAL_vertex[n].neigh;

            i = GLOBAL_vertex[n].px;
            j = GLOBAL_vertex[n].py;
            k = GLOBAL_vertex[n].pz;

            if (neighbor & 1)
            {
                m++;
                l   = GLOBAL_atom_index[IndexVect(i - 1, j, k)];
                nx += GLOBAL_vertex[l].x;
                ny += GLOBAL_vertex[l].y;
                nz += GLOBAL_vertex[l].z;
            }

            if (neighbor & 2)
            {
                m++;
                l   = GLOBAL_atom_index[IndexVect(i + 1, j, k)];
                nx += GLOBAL_vertex[l].x;
                ny += GLOBAL_vertex[l].y;
                nz += GLOBAL_vertex[l].z;
            }

            if (neighbor & 4)
            {
                m++;
                l   = GLOBAL_atom_index[IndexVect(i, j - 1, k)];
                nx += GLOBAL_vertex[l].x;
                ny += GLOBAL_vertex[l].y;
                nz += GLOBAL_vertex[l].z;
            }

            if (neighbor & 8)
            {
                m++;
                l   = GLOBAL_atom_index[IndexVect(i, j + 1, k)];
                nx += GLOBAL_vertex[l].x;
                ny += GLOBAL_vertex[l].y;
                nz += GLOBAL_vertex[l].z;
            }

            if (neighbor & 16)
            {
                m++;
                l   = GLOBAL_atom_index[IndexVect(i, j, k - 1)];
                nx += GLOBAL_vertex[l].x;
                ny += GLOBAL_vertex[l].y;
                nz += GLOBAL_vertex[l].z;
            }

            if (neighbor & 32)
            {
                m++;
                l   = GLOBAL_atom_index[IndexVect(i, j, k + 1)];
                nx += GLOBAL_vertex[l].x;
                ny += GLOBAL_vertex[l].y;
                nz += GLOBAL_vertex[l].z;
            }

            // update the position
            GLOBAL_vertex[n].x = nx / (float)m;
            GLOBAL_vertex[n].y = ny / (float)m;
            GLOBAL_vertex[n].z = nz / (float)m;
        }
    }
    finish = clock();

    // printf("   Smooth the quad meshes: CPU Time = %f seconds
    // \n\n",(double)(finish-begin)/CLOCKS_PER_SEC);

    // Allocate memory
    std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh);

    // write vertices
    for (int i = 0; i < GLOBAL_vert_num; i++)
    {
        float x = GLOBAL_vertex[i].x * span[0] + orig[0];
        float y = GLOBAL_vertex[i].y * span[1] + orig[1];
        float z = GLOBAL_vertex[i].z * span[2] + orig[2];

        mesh->insert<1>({i}, SMVertex({x, y, z}));
    }

    // write triangles
    float angle, angle1, angle2;

    for (i = 0; i < GLOBAL_quad_num; i++)
    {
        a = GLOBAL_quads[i].a;
        b = GLOBAL_quads[i].b;
        c = GLOBAL_quads[i].c;
        d = GLOBAL_quads[i].d;

        angle1 = -999.0;
        angle2 = -999.0;
        angle  = GetAngle(a, b, c);

        if (angle > angle1)
        {
            angle1 = angle;
        }
        angle = GetAngle(a, d, c);

        if (angle > angle1)
        {
            angle1 = angle;
        }
        angle = GetAngle(c, a, b);

        if (angle > angle1)
        {
            angle1 = angle;
        }
        angle = GetAngle(c, a, d);

        if (angle > angle1)
        {
            angle1 = angle;
        }

        angle = GetAngle(b, a, d);

        if (angle > angle2)
        {
            angle2 = angle;
        }
        angle = GetAngle(b, c, d);

        if (angle > angle2)
        {
            angle2 = angle;
        }
        angle = GetAngle(d, a, b);

        if (angle > angle2)
        {
            angle2 = angle;
        }
        angle = GetAngle(d, b, c);

        if (angle > angle2)
        {
            angle2 = angle;
        }

        if (angle1 <= angle2)
        {
            mesh->insert<3>({a, b, c});
            mesh->insert<3>({a, c, d});
        }
        else
        {
            mesh->insert<3>({a, b, d});
            mesh->insert<3>({b, c, d});
        }
    }

    /* write to disk
       FILE *fout;
       if ((fout=fopen("output1.off", "wb"))==NULL){
       printf("write error...\n");
       exit(0);
       };

       fprintf(fout, "OFF\n");
       fprintf(fout, "%d %d
          %d\n",GLOBAL_vert_num,quad_num*2,GLOBAL_vert_num+quad_num*2-2);
       for (i = 0; i < GLOBAL_vert_num; i++)
       fprintf(fout, "%f %f %f
          \n",surfmesh->vertex[i].x,surfmesh->vertex[i].y,surfmesh->vertex[i].z);

       for (i = 0; i < quad_num*2; i++) {
       fprintf(fout, "3 %d %d
          %d\n",surfmesh->face[i].a,surfmesh->face[i].b,surfmesh->face[i].c);
       }
       fclose(fout);
     */

    free(AllSeeds);
    free(GLOBAL_segment_index);
    free(GLOBAL_atom_index);
    free(min_heap->x);
    free(min_heap->y);
    free(min_heap->z);
    free(min_heap->seed);
    free(min_heap->dist);
    free(min_heap);
    free(GLOBAL_vertex);
    free(GLOBAL_quads);

    compute_orientation(*mesh);
    return mesh;
}

float GetAngle(int a, int b, int c)
{
    float ax, ay, az;
    float bx, by, bz;
    float dist;

    ax   = GLOBAL_vertex[b].x - GLOBAL_vertex[a].x;
    ay   = GLOBAL_vertex[b].y - GLOBAL_vertex[a].y;
    az   = GLOBAL_vertex[b].z - GLOBAL_vertex[a].z;
    dist = sqrt(ax * ax + ay * ay + az * az);

    if (dist > 0)
    {
        ax /= dist;
        ay /= dist;
        az /= dist;
    }
    bx   = GLOBAL_vertex[c].x - GLOBAL_vertex[a].x;
    by   = GLOBAL_vertex[c].y - GLOBAL_vertex[a].y;
    bz   = GLOBAL_vertex[c].z - GLOBAL_vertex[a].z;
    dist = sqrt(bx * bx + by * by + bz * bz);

    if (dist > 0)
    {
        bx /= dist;
        by /= dist;
        bz /= dist;
    }

    return ax * bx + ay * by + az * bz;
}

char CheckManifold(int i, int j, int k)
{
    char manifold, nonmanifold;
    int  m, n, l;

    manifold = 1;

    for (l = k - 1; l <= k + 1; l++)
    {
        if (!manifold)
        {
            break;
        }

        for (n = j - 1; n <= j + 1; n++)
        {
            if (!manifold)
            {
                break;
            }

            for (m = i - 1; m <= i + 1; m++)
            {
                if ((m != i) || (n != j) || (l != k))
                {
                    if (GLOBAL_segment_index[IndexVect(m, n, l)] == MaxVal)
                    {
                        nonmanifold = 1;

                        if ((m != i) && (GLOBAL_segment_index[IndexVect(m, j, k)] == MaxVal))
                        {
                            nonmanifold = 0;
                        }

                        if ((n != j) && (GLOBAL_segment_index[IndexVect(i, n, k)] == MaxVal))
                        {
                            nonmanifold = 0;
                        }

                        if ((l != k) && (GLOBAL_segment_index[IndexVect(i, j, l)] == MaxVal))
                        {
                            nonmanifold = 0;
                        }

                        if (nonmanifold)
                        {
                            manifold = 0;
                        }
                    }
                }

                if (!manifold)
                {
                    break;
                }
            }
        }
    }

    return manifold;
}

int CheckFaceCorner(float x, float y, float z)
{
    int m, n, l;
    int a;

    m = (int)x;
    n = (int)y;
    l = (int)z;

    if (GLOBAL_atom_index[IndexVect(m, n, l)] < 0)
    {
        GLOBAL_vertex[GLOBAL_vert_num].x  = x;
        GLOBAL_vertex[GLOBAL_vert_num].y  = y;
        GLOBAL_vertex[GLOBAL_vert_num].z  = z;
        GLOBAL_vertex[GLOBAL_vert_num].px = m;
        GLOBAL_vertex[GLOBAL_vert_num].py = n;
        GLOBAL_vertex[GLOBAL_vert_num].pz = l;
        GLOBAL_atom_index[IndexVect(m, n, l)] = GLOBAL_vert_num;
        a = GLOBAL_vert_num;
        GLOBAL_vert_num++;
    }
    else
    {
        a = GLOBAL_atom_index[IndexVect(m, n, l)];
    }

    return a;
}

FLT2VECT FindIntersection(int n, int m, int j, int k, ATOM *atom_list)
{
    FLT2VECT intersect;
    int      i;
    char     reorder;
    double   cx1, cy1, cz1;
    double   cx2, cy2, cz2;
    double   dist, radius1, radius2;
    double   cos_alpha;
    double   ax, ay, az;
    double   cx, cy, cz;


    n--;
    m--;
    reorder = 0;

    // always make the first atom on the left
    // and the second on the right
    if (atom_list[n].x > atom_list[m].x)
    {
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

    if (cx1 == cx2)
    {
        if ((radius1 * radius1 - (k - cz1) * (k - cz1) - (j - cy1) * (j - cy1)) >=
            (radius2 * radius2 - (k - cz2) * (k - cz2) - (j - cy2) * (j - cy2)))
        {
            intersect.x = (float)MaxVal;
            intersect.y = -(float)MaxVal;
        }
        else
        {
            intersect.x = -(float)MaxVal;
            intersect.y = (float)MaxVal;
        }
    }
    else
    {
        dist = sqrt((cx2 - cx1) * (cx2 - cx1) + (cy2 - cy1) * (cy2 - cy1) + (cz2 - cz1) * (cz2 - cz1));
        cos_alpha = (radius1 * radius1 + dist * dist - radius2 * radius2) / (2.0 * radius1 * dist);
        ax   = (cx2 - cx1) / dist;
        ay   = (cy2 - cy1) / dist;
        az   = (cz2 - cz1) / dist;
        dist = radius1 * cos_alpha;
        cx   = dist * ax + cx1;
        cy   = dist * ay + cy1;
        cz   = dist * az + cz1;
        dist = cx * ax + (cy - j) * ay + (cz - k) * az;

        if (reorder == 0)
        {
            intersect.x = (float)(dist / ax);
            intersect.y = (float)MaxVal;
        }
        else
        {
            intersect.x = -(float)MaxVal;
            intersect.y = (float)(dist / ax);
        }
    }
    return intersect;
}

int ExtractSAS(int atom_num, ATOM *atom_list)
{
    int      i, j, k;
    int      m, n, l;
    int      dim[3], c[3];
    int      amax[3], amin[3];
    float    radius;
    float    x, y, z;
    FLT2VECT intersect;


    dim[0] = GLOBAL_xdim;
    dim[1] = GLOBAL_ydim;
    dim[2] = GLOBAL_zdim;


    for (m = 0; m < atom_num; m++)
    {
        radius = atom_list[m].radius;

        // compute the dataset coordinates of the atom's center
        c[0] = (int)(atom_list[m].x + 0.5);
        c[1] = (int)(atom_list[m].y + 0.5);
        c[2] = (int)(atom_list[m].z + 0.5);

        // then compute the bounding box of the atom
        for (j = 0; j < 3; j++)
        {
            int tmp;
            tmp = (int)(c[j] - radius - 1);
            tmp = (tmp < 0) ? 0 : tmp;
            amin[j] = tmp;
            tmp = (int)(c[j] + radius + 1);
            tmp = (tmp > (dim[j] - 1)) ? (dim[j] - 1) : tmp;
            amax[j] = tmp;
        }

        // begin blurring kernel
        radius = radius * radius;

        for (k = amin[2]; k <= amax[2]; k++)
        {
            for (j = amin[1]; j <= amax[1]; j++)
            {
                for (i = amin[0]; i <= amax[0]; i++)
                {
                    x = i;
                    y = j;
                    z = k;

                    if ((x - atom_list[m].x) * (x - atom_list[m].x) +
                        (y - atom_list[m].y) * (y - atom_list[m].y) +
                        (z - atom_list[m].z) * (z - atom_list[m].z) <= radius)
                    {
                        if (GLOBAL_atom_index[IndexVect(i, j, k)] > 0)
                        {
                            intersect = FindIntersection(GLOBAL_atom_index[IndexVect(i, j, k)], m + 1, j, k, atom_list);

                            if ((i >= intersect.x) && (i <= intersect.y))
                            {
                                GLOBAL_atom_index[IndexVect(i, j, k)] = m + 1;
                            }
                        }
                        else
                        {
                            GLOBAL_atom_index[IndexVect(i, j, k)] = m + 1;
                        }
                    }
                }
            }
        }
    }


    // find voxels on the border
    int total = 0;

    for (l = 1; l < GLOBAL_zdim - 1; l++)
    {
        for (n = 1; n < GLOBAL_ydim - 1; n++)
        {
            for (m = 1; m < GLOBAL_xdim - 1; m++)
            {
                if (GLOBAL_atom_index[IndexVect(m, n, l)])
                {
                    int count = 0;

                    // Look at neighbor voxels
                    for (k = std::max(l - 1, 0); k <= std::min(l + 1, GLOBAL_zdim - 1); k++)
                    {
                        for (j = std::max(n - 1, 0); j <= std::min(n + 1, GLOBAL_ydim - 1); j++)
                        {
                            for (i = std::max(m - 1, 0); i <= std::min(m + 1, GLOBAL_xdim - 1); i++)
                            {
                                if ((((i == m) && (j == n)) || ((i == m) && (k == l)) || ((k == l) && (j == n))) &&
                                    (GLOBAL_atom_index[IndexVect(i, j, k)] == 0))
                                {
                                    count = 1;
                                }
                            }
                        }
                    }

                    if (count)
                    {
                        GLOBAL_atom_index[IndexVect(m, n, l)] = -GLOBAL_atom_index[IndexVect(m, n, l)];
                        total++;
                    }
                }
            }
        }
    }


    return total;
}

} // end namespace gamer
