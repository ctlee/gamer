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

#ifndef TETSMOOTH_H_UWM
#define TETSMOOTH_H_UWM

#include "DataStruct.h"

// declaration of functions
int			Usage(void);
int			DisplayMeshInfo(void);

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

void		PerformAllSmooth(void);
int			PerformSingleStep(const int verIndex);

DBL3DVECT	InnerODT_SingleTet(const DBL3DVECT& curV,const DBL3DVECT& NB1,const DBL3DVECT& NB2,const DBL3DVECT& NB3);
int			Smooth_3DOriginalODT(const int verIndex);
int			SmoothBdyVer(const int verIndex);
EIGENVECT	GetEigenVector(verElement* curVer,double eigenValues[3]);
int			tmpGet_EFG(const int verIndex,double& E,double& F,double& G,double& H,double& I,DBL3DVECT& dirS,DBL3DVECT& dirT);	
int			Smooth_2DODT_VolFixedOptimization(const int curVer);
int			getBDY_ODT_FeaturePres_Paras(const int verIndex,double& A,double& B,const DBL3DVECT& dirMainCurv);
int			Smooth_BDYODT_FeaturePreserving(const int verIndex);

void		ComputeFaceNormal( const verElement* const Ver0,\
							  const verElement* const Ver1,\
							  const verElement* const Ver2,\
							  const verElement* const Ver3,\
							  double faceNormal[3],bool PointInside = true);
EIGENVECT	EigenVectors_Jocobi(double Matrix[3][3],double eigenValues[3]);

bool		CheckValid_byPosVol(const int verIndex);
double		Area(const DBL3DVECT& P1,const DBL3DVECT& P2,const DBL3DVECT& P3);
double		Area(const verElement* v1, const verElement* v2,const verElement* v3);
double		Volume(const DBL3DVECT& P1,const DBL3DVECT& P2,const DBL3DVECT& P3,const DBL3DVECT& P4);
double		SignedVolume(const DBL3DVECT& P1,const DBL3DVECT& P2,const DBL3DVECT& P3,const DBL3DVECT& P4);
double		SignedVolume(const tetElement* curTet);
void		InitVolumeSigns(void);

int			ComputeDihedralAngles_All(void);// to statistic the distribution 
double		ComputeLocalAngle(const int verIndex,double& localMaxAngle);
double		ComputeMinDiAngle(const DBL3DVECT& P1,const DBL3DVECT& P2,const DBL3DVECT& P3,const DBL3DVECT& P4);
int			ComputeDiAngle_SigleTet(const tetElement* curTet, double& minAngle,double& maxAngle);
double		ComputeDiAngle(const verElement* P,const verElement* Q,const verElement* r1,const verElement* r2);
double		ComputeDiAngle(const DBL3DVECT& P,const DBL3DVECT& Q,const DBL3DVECT& r1,const DBL3DVECT& r2);

DBL3DVECT	VerElement2DBLVECT(const verElement* ver);
void		DBLVECT2VerElement(const DBL3DVECT &DBL, verElement* &ver);

#endif