/*****************************************************************************
 * Tetrahedral Mesh Smoothing
 * Copyright (C) 2011--  Zhanheng Gao and Zeyun Yu
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 675 Mass Ave, Cambridge, MA 02139, USA.
 * ***************************************************************************/

#ifndef TET2TET3_H_UWM
#define TET2TET3_H_UWM

#include "DataStruct.h"

double		Area(const DBL3DVECT& P1,const DBL3DVECT& P2,const DBL3DVECT& P3);
DBL3DVECT	VerElement2DBLVECT(const verElement* ver);
void		DBLVECT2VerElement(const DBL3DVECT &DBL, verElement* &ver);

int			InitBasicData(void);
int			LoadData(const char* filename);
int			SaveData(const char* filename);
int			ComputeTetInfo_forVer(void);
int			ComputeBdyInfo_forTet(void);
int			CorrectTetInfo(void);
DBL3DVECT	ComputeFaceNormal(const triElement* face,bool PointInside = true);
void		ComputeAllFaceNormal(void);
int			ComputeBdyVerNormal(void);
int			ReleaseAll(void);

int			reComputeBdyInfoForTet(void);
bool		SearchFaceIn_AllTets(const int verIndex1,const int verIndex2,const int verIndex3,const int skipTetIndex);
int			GenerateTetFaces(void);
int			SearchtetFACES(int curTetIndex,int ivTri[3],int ivOther,std::vector<triElement> &faceVector);
int			ComputeTriNormal3D(triElement& tri);

#endif