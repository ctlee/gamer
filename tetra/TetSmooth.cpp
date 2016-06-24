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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "TetSmooth.h"

// global variables
int							m_numVer = 0;
std::vector<verElement*>	verList;

int							m_numTetFACES = 0;
std::vector<triElement>		tetFACES;

int							m_numTet = 0;
std::vector<tetElement>		tetMESH;
std::vector<int>			m_VolumeSignes;

double						m_minMinAngle = 181;
double						m_maxMaxAngle = -1;

double						ori_minX = 0,ori_maxX = 0;
double						ori_minY = 0,ori_maxY = 0;
double						ori_minZ = 0,ori_maxZ = 0;
double						m_oriMetric = 0;

int Usage(void)
{
	printf("Usage: TetSmooth <input Tet3 file> <output Tet3 file> [iteration times (1 for default)]\n");
	return 0;
}

int DisplayMeshInfo(void)
{
	printf("The min and max dihedral angle in current mesh are: %.2f and %.2f (degree).\n",m_minMinAngle,m_maxMaxAngle);
	return 0;
}
int main(int argc,char** argv)
{
	using namespace std;
	if (argc < 3)
	{
		Usage();
		return 0;
	}

	const char* inputFileName = argv[1];
	const char* outputFileName = argv[2];
	float numIterations = 0;
	if (argc >= 4)
		numIterations = atof(argv[3]);
	else
		numIterations = 1;

	printf("Initializing basic data....\n");
	InitBasicData();

	printf("Loading Data....\n");
	LoadData(inputFileName);
	DisplayMeshInfo();

	printf("Smoothing....\n");
	int k = 0;
	for (k = 0; k < numIterations; k ++)
	{
		PerformAllSmooth();
		ComputeDihedralAngles_All();
		printf("Finished smoothing %d times, current min and max angles are %.2f and %.2f (degree).\n",k + 1,m_minMinAngle,m_maxMaxAngle);
	}

	printf("Smoothing finished.\n");
	ComputeDihedralAngles_All();
	DisplayMeshInfo();

	printf("Saving result data....\n");
	SaveData(outputFileName);

	printf("Releasing memory....\n");
	ReleaseAll();

	return 0;
}

void PerformAllSmooth(void)
{
	
	if (verList.size() == 0 || tetFACES.size() == 0)
	{
		return;
	}	

	int numVers = m_numVer;

	ComputeDihedralAngles_All();
	const double minAngle_GlobalMesh = m_minMinAngle;
	const double maxAngle_GlobalMesh = m_maxMaxAngle;
	double localMinAngle,newLocalMinAngle;
	double localMaxAngle,newLocalMaxAngle;

//	DWORD thisClock = GetTickCount();
	DBL3DVECT OriCoord;	
	int iver = 0;
	for ( iver = 0; iver < numVers; iver ++)
	{
		localMinAngle = ComputeLocalAngle(iver,localMaxAngle);

		OriCoord = VerElement2DBLVECT(verList[iver]);
		if (0 != PerformSingleStep(iver))// the coordinate of verList[iver] will be modified in this function.
		{
			DBLVECT2VerElement(OriCoord,verList[iver]);
			continue;
		}

		if (!CheckValid_byPosVol(iver))// the new position of current node will lead to invalid elements according to the positive volume rule
		{
			{
				DBLVECT2VerElement(OriCoord,verList[iver]);
				continue;
			}
		}		

		newLocalMinAngle = ComputeLocalAngle(iver,newLocalMaxAngle);
		if (newLocalMinAngle < localMinAngle)
		{// guarantee the improvement of the min angle 
			DBLVECT2VerElement(OriCoord,verList[iver]);
		}
	}	

//	thisClock = GetTickCount() - thisClock;
}

int PerformSingleStep(const int verIndex)
{	
	verElement* curVer = verList[verIndex];
	
	if (curVer->m_type == PTTYPE_BDY )
	{
		return SmoothBdyVer(verIndex);		
	}
	else
	{
		return Smooth_3DOriginalODT(verIndex);// smooth in ODT manner for inner nodes		
	}
	
	return 0;
}

DBL3DVECT InnerODT_SingleTet(const DBL3DVECT& curV,const DBL3DVECT& NB1,const DBL3DVECT& NB2,const DBL3DVECT& NB3)
{
	DBL3DVECT SingleTerm;
	
	DBL3DVECT faceNormal;
	DBL3DVECT RefVec;
	double faceArea;
	double edgeLengthSquare[3];
	double totalLengthSquare;
	
	RefVec.x = curV.x - NB1.x;
	RefVec.y = curV.y - NB1.y;
	RefVec.z = curV.z - NB1.z;
	
	faceNormal = (NB2 - NB1) / (NB3 - NB1);
	faceArea = faceNormal.Len() / 2;
	faceNormal /= faceNormal.Len();

	if ((faceNormal * RefVec) < 0)
	{
		faceNormal = -faceNormal;
	}	
	
	edgeLengthSquare[0] = pow((curV - NB1).Len(),2);
	edgeLengthSquare[1] = pow((curV - NB2).Len(),2);
	edgeLengthSquare[2] = pow((curV - NB3).Len(),2);
	totalLengthSquare = edgeLengthSquare[0] + edgeLengthSquare[1] + edgeLengthSquare[2];
	
	SingleTerm.x = faceArea * totalLengthSquare * faceNormal.x / 3;
	SingleTerm.y = faceArea * totalLengthSquare * faceNormal.y / 3;
	SingleTerm.z = faceArea * totalLengthSquare * faceNormal.z / 3;	
	
	return SingleTerm;
}

int Smooth_3DOriginalODT(const int verIndex)
{
	verElement* curVer = verList[verIndex];
	
	int numAdjTet = curVer->m_numTet;
	int *adjTet = curVer->m_tets;	
	
	DBL3DVECT Q[4];	
	Q[0] = VerElement2DBLVECT(curVer);
	
	DBL3DVECT V_star;	
	double tetVolume;
	double totalVolume = 0;
	tetElement* curTet;
	
	int indexOtherVer[3];
	
	DBL3DVECT  SingleTerm;
	V_star.x = V_star.y = V_star.z = 0;

	for (int itet = 0; itet < numAdjTet; itet ++)
	{
		curTet = &(tetMESH[adjTet[itet]]);// adjTet[itet] is just the exist tet.
		indexOtherVer[0] = curVer->m_subscriptIn_tet_2nd[itet];
		indexOtherVer[1] = curVer->m_subscriptIn_tet_3rd[itet];
		indexOtherVer[2] = curVer->m_subscriptIn_tet_4th[itet];		
		indexOtherVer[0] = curTet->m_v[indexOtherVer[0]];
		indexOtherVer[1] = curTet->m_v[indexOtherVer[1]];
		indexOtherVer[2] = curTet->m_v[indexOtherVer[2]];				
		Q[1] = VerElement2DBLVECT(verList[indexOtherVer[0]]);
		Q[2] = VerElement2DBLVECT(verList[indexOtherVer[1]]);
		Q[3] = VerElement2DBLVECT(verList[indexOtherVer[2]]);		
		
		tetVolume = Volume(Q[0],Q[1],Q[2],Q[3]);		
		totalVolume += tetVolume;		
		
		SingleTerm = InnerODT_SingleTet(Q[0],Q[1],Q[2],Q[3]);
		
		V_star.x += SingleTerm.x;
		V_star.y += SingleTerm.y;
		V_star.z += SingleTerm.z;
	}
	
	curVer->m_x -= V_star.x / 2 / totalVolume;
	curVer->m_y -= V_star.y / 2 / totalVolume;
	curVer->m_z -= V_star.z / 2 / totalVolume;
	
	return 0;
}

int SmoothBdyVer(const int verIndex)
{
	verElement* curVer = verList[verIndex];
	
	double ka = 0.1;
	double eigenValues[3];
	GetEigenVector(curVer,eigenValues);
	if (fabs(eigenValues[2] / eigenValues[0]) > ka)
		;
	else if (fabs(eigenValues[1] / eigenValues[0]) > ka && fabs(eigenValues[2] / eigenValues[0]) < ka)// edge
		return Smooth_BDYODT_FeaturePreserving(verIndex);
	else
		return Smooth_2DODT_VolFixedOptimization(verIndex);
	
	return 0;
}

int Smooth_2DODT_VolFixedOptimization(const int verIndex)
{
	verElement* curVer = verList[verIndex];
	
	double E,F,G,H,I;
	DBL3DVECT dirS,dirT;
	if (0 != tmpGet_EFG(verIndex,E,F,G,H,I,dirS,dirT))
	{
		return -1;
	}
	
	double Lower = 4 * E * F - G * G;
	double Up1 = G * I - 2 * F * H;
	double Up2 = G * H - 2 * E * I;
	
	double coordS = Up1 / Lower;
	double coordT = Up2 / Lower;
	
	
	curVer->m_x += coordS * dirS.x + coordT * dirT.x;
	curVer->m_y += coordS * dirS.y + coordT * dirT.y;
	curVer->m_z += coordS * dirS.z + coordT * dirT.z;
	
	return 0;
}

int tmpGet_EFG(const int verIndex,double& E,double& F,double& G,double& H,double& I,DBL3DVECT& coorS,DBL3DVECT& coorT)
{
	verElement* curVer = verList[verIndex];

	DBL3DVECT curVDBL;
	double triArea;
	DBL3DVECT triNormal;
	double norMag;
	DBL3DVECT normalVector;

	int numTet = curVer->m_numTet;
	int* adjTet = curVer->m_tets;
	tetElement* curTet;

	int index_NbrsIn_Tet[3];
	DBL3DVECT nbrVertex[3];// the three neighbors of curV in every tet
	DBL3DVECT nbrDifVector[3];// the three difference vectors from curVer to its three tet neighbors
	double volOmega_star;
	curVDBL = VerElement2DBLVECT(curVer);

	triElement* curTri;
	int index_NbrsIn_Tri[2];
	DBL3DVECT	nbrX[2];// the two DIFFERENCE vectors of from curV to its boundary adjacent neighbors. conunter-clockwise 
	DBL3DVECT   oppoNbr;// the opposite vertex of current vertex in a tet which has a face consisted of current vertex and its two neighbors
	DBL3DVECT   sum_NbrX;
	DBL3DVECT   cross_NbrX;

	// first, compute the volume of Omega_0 (which is also Omega_star)
	// and get the normal vector which determinans the tangent plane
	volOmega_star = 0;
	normalVector.x = normalVector.y = normalVector.z = 0;

	DBL3DVECT tetFaceNormals[4];
	DBL3DVECT tetFaceNormalSum(0,0,0);
	int itet = 0;
	for (itet = 0; itet < numTet; itet ++)
	{
		curTet = &(tetMESH[curVer->m_tets[itet]]);
		index_NbrsIn_Tet[0] = curVer->m_subscriptIn_tet_2nd[itet];
		index_NbrsIn_Tet[1] = curVer->m_subscriptIn_tet_3rd[itet];
		index_NbrsIn_Tet[2] = curVer->m_subscriptIn_tet_4th[itet];		
		index_NbrsIn_Tet[0] = curTet->m_v[index_NbrsIn_Tet[0]];
		index_NbrsIn_Tet[1] = curTet->m_v[index_NbrsIn_Tet[1]];
		index_NbrsIn_Tet[2] = curTet->m_v[index_NbrsIn_Tet[2]];		

		nbrVertex[0] = VerElement2DBLVECT(verList[index_NbrsIn_Tet[0]]);
		nbrVertex[1] = VerElement2DBLVECT(verList[index_NbrsIn_Tet[1]]);
		nbrVertex[2] = VerElement2DBLVECT(verList[index_NbrsIn_Tet[2]]);		
		volOmega_star += Volume(curVDBL,nbrVertex[0],nbrVertex[1],nbrVertex[2]);

		nbrDifVector[0] = nbrVertex[0] - curVDBL;
		nbrDifVector[0] /= nbrDifVector[0].Len();
		triArea = Area(nbrVertex[0],nbrVertex[1],nbrVertex[2]);
		triNormal = (nbrVertex[1] - nbrVertex[0]) / (nbrVertex[2] - nbrVertex[0]);
		norMag = triNormal.Len();
		triNormal = triNormal / norMag;		
		if (triNormal * nbrDifVector[0] > 1e-10)
		{
			triNormal = -triNormal;
		}
		normalVector = normalVector +  triNormal * triArea;	
	}

	const double NormalLen = norMag = normalVector.Len();
	normalVector = normalVector / norMag;// generate the unit normal vector

	// second, get the s and t coordinate system directional vectors
	if (fabs(normalVector.x) > 1e-10)
	{
		coorS.y = 1;
		coorS.z = 0;
		coorS.x = - normalVector.y / normalVector.x;
	}
	else
	{
		coorS.x = 1;
		coorS.z = 0;
		coorS.y = 0;
	}	

	norMag = coorS.Len();
	coorS = coorS / norMag;
	
	coorT = normalVector / coorS;
	norMag = coorT.Len();
	coorT = coorT / norMag;	
	if (fabs(normalVector * coorS > 1e-10 || fabs(normalVector * coorT) > 1e-10))
	{
		perror("Computing CoorS and CoorT error.\n");
	}
	normalVector = normalVector * NormalLen;// rescale the unit normal vector to the original length

	// third, compute the E,F,G
	E = F = volOmega_star / 4;
	G = 0;
	DBL3DVECT refVect;
	DBL3DVECT chkNormal(0,0,0);
	DBL3DVECT chkCrossPdt;
	double    chkDiff;
	int numOmittedFaces = 0;
	int numValidFaces = 0;	

	int idxCurTet;
	tetFaceNormalSum.x = tetFaceNormalSum.y = tetFaceNormalSum.z = 0;
 	for (itet = 0; itet < numTet; itet ++)
	{
		idxCurTet = curVer->m_tets[itet];
		curTet = &(tetMESH[idxCurTet]);
		for (int iFace = 0; iFace < 4; iFace ++)
		{
			curTri = &(tetFACES[curTet->m_faces[iFace]]);

			if (verIndex != curTri->m_v[0] && verIndex != curTri->m_v[1] && verIndex != curTri->m_v[2])
				continue;
			if (curTri->m_type != FACETYPE_BDY)
			{
				numOmittedFaces ++;
				continue;
			}
			if (curTri->m_tet[1] != -1 || curTri->m_ssb[1] != -1)
			{
				perror("Boundary face error.");
			}
			
			numValidFaces ++;
			if (verIndex == curTri->m_v[0])
			{
				index_NbrsIn_Tri[0] = curTri->m_v[1];
				index_NbrsIn_Tri[1] = curTri->m_v[2];
			}
			else if (verIndex == curTri->m_v[1])
			{
				index_NbrsIn_Tri[0] = curTri->m_v[0];
				index_NbrsIn_Tri[1] = curTri->m_v[2];
			}
			else if (verIndex == curTri->m_v[2])
			{
				index_NbrsIn_Tri[0] = curTri->m_v[0];
				index_NbrsIn_Tri[1] = curTri->m_v[1];
			}
			else
			{
				perror("Error, \n");
			}		

			nbrX[0] = VerElement2DBLVECT(verList[index_NbrsIn_Tri[0]]);
			nbrX[1] = VerElement2DBLVECT(verList[index_NbrsIn_Tri[1]]);

			triArea = Area(curVDBL,nbrX[0],nbrX[1]);

			nbrX[0] = nbrX[0] - curVDBL;		
			nbrX[1] = nbrX[1] - curVDBL;		
			sum_NbrX = nbrX[0] + nbrX[1];		
			cross_NbrX = nbrX[0] / nbrX[1];
			if (curTri->m_tet[0] == idxCurTet)
			{
				oppoNbr = VerElement2DBLVECT(verList[curTet->m_v[iFace]]);
			}
			else
			{
				perror("Error, face tet index.\n");
			}

			refVect = oppoNbr - curVDBL;
			refVect /= refVect.Len();	
 			triNormal = cross_NbrX / cross_NbrX.Len();
			if (triNormal * refVect > 1e-10)
			{
				cross_NbrX = -cross_NbrX;

				triNormal = -triNormal;
			}
			
			chkNormal = chkNormal + triNormal * triArea;

			E -= (coorS * sum_NbrX) * (coorS * cross_NbrX) / 60;
			F -= (coorT * sum_NbrX) * (coorT * cross_NbrX) / 60;
			G -= ((coorS * sum_NbrX) * (coorT * cross_NbrX) + (coorT * sum_NbrX) * (coorS * cross_NbrX));
		}
	}
	G /= 60;

	chkCrossPdt = chkNormal - normalVector;
	chkDiff = chkCrossPdt.Len();
	if (chkDiff > 1e-10)// means there is some self-intersections occurs. i.e. this can be looked as a way of checking inverted elements. May 23, 2011
	{
		return -1;
	}

	if (E > 0 && 4 * E * F - G * G > 0)
	{
	}
	else
	{
 		perror("There is non Possitive Definite Matrix!\n");
		return -1;
	}	

	// fourth, compute H, and I
	H = I = 0;
	for (itet = 0; itet < numTet; itet ++)
	{
		curTet = &(tetMESH[curVer->m_tets[itet]]);
		index_NbrsIn_Tet[0] = curVer->m_subscriptIn_tet_2nd[itet];
		index_NbrsIn_Tet[1] = curVer->m_subscriptIn_tet_3rd[itet];
		index_NbrsIn_Tet[2] = curVer->m_subscriptIn_tet_4th[itet];		
		index_NbrsIn_Tet[0] = curTet->m_v[index_NbrsIn_Tet[0]];
		index_NbrsIn_Tet[1] = curTet->m_v[index_NbrsIn_Tet[1]];
		index_NbrsIn_Tet[2] = curTet->m_v[index_NbrsIn_Tet[2]];		
		
		nbrVertex[0] = VerElement2DBLVECT(verList[index_NbrsIn_Tet[0]]);
		nbrVertex[1] = VerElement2DBLVECT(verList[index_NbrsIn_Tet[1]]);
		nbrVertex[2] = VerElement2DBLVECT(verList[index_NbrsIn_Tet[2]]);		


		nbrDifVector[0] = nbrVertex[0] - curVDBL;
		nbrDifVector[0] /= nbrDifVector[0].Len();
		triArea = Area(nbrVertex[0],nbrVertex[1],nbrVertex[2]);
		triNormal = (nbrVertex[1] - nbrVertex[0]) / (nbrVertex[2] - nbrVertex[0]);
		norMag = sqrt(triNormal.x * triNormal.x + triNormal.y * triNormal.y + triNormal.z * triNormal.z);
		triNormal = triNormal / norMag;		
		if (triNormal * nbrDifVector[0] > 1e-10)
		{
			triNormal = -triNormal;
		}
		H += triArea * (coorS * triNormal) * (pow((nbrVertex[0] - curVDBL).Len(),2) + \
											  pow((nbrVertex[1] - curVDBL).Len(),2) + 
											  pow((nbrVertex[2] - curVDBL).Len(),2));
		I += triArea * (coorT * triNormal) * (pow((nbrVertex[0] - curVDBL).Len(),2) + \
											  pow((nbrVertex[1] - curVDBL).Len(),2) + \
											  pow((nbrVertex[2] - curVDBL).Len(),2));
	}

	H /= 12;
	I /= 12;

	for (itet = 0; itet < numTet; itet ++)
	{
		idxCurTet = curVer->m_tets[itet];
		curTet = &(tetMESH[idxCurTet]);
		for (int iFace = 0; iFace < 4; iFace ++)
		{
			curTri = &(tetFACES[curTet->m_faces[iFace]]);

			if (curTri->m_type != FACETYPE_BDY)
				continue;
			if (verIndex != curTri->m_v[0] && verIndex != curTri->m_v[1] && verIndex != curTri->m_v[2])
				continue;

			if (verIndex == curTri->m_v[0])
			{
				index_NbrsIn_Tri[0] = curTri->m_v[1];
				index_NbrsIn_Tri[1] = curTri->m_v[2];
			}
			else if (verIndex == curTri->m_v[1])
			{
				index_NbrsIn_Tri[0] = curTri->m_v[0];
				index_NbrsIn_Tri[1] = curTri->m_v[2];
			}
			else if (verIndex == curTri->m_v[2])
			{
				index_NbrsIn_Tri[0] = curTri->m_v[0];
				index_NbrsIn_Tri[1] = curTri->m_v[1];
			}
			else
			{
				perror("Error.");
			}
			
			nbrX[0] = VerElement2DBLVECT(verList[index_NbrsIn_Tri[0]]);
			nbrX[1] = VerElement2DBLVECT(verList[index_NbrsIn_Tri[1]]);
			nbrX[0] = nbrX[0] - curVDBL;		
			nbrX[1] = nbrX[1] - curVDBL;		
			cross_NbrX = nbrX[0] / nbrX[1];
			
			if (idxCurTet == curTri->m_tet[0])
			{
				oppoNbr = VerElement2DBLVECT(verList[curTet->m_v[iFace]]);
			}
			else
			{
				perror("Error, face tet index.\n");
			}

			oppoNbr = oppoNbr - curVDBL;	
			oppoNbr /= oppoNbr.Len();
			triNormal = cross_NbrX / cross_NbrX.Len();
			if (triNormal * oppoNbr > 1e-10)
			{
				cross_NbrX = -cross_NbrX;
			}

			double tmp = pow(nbrX[0].Len(),2);

			H -= (pow(nbrX[0].Len(),2) + pow(nbrX[1].Len(),2) + (nbrX[0] * nbrX[1])) * (coorS * cross_NbrX) / 60;
			I -= (pow(nbrX[0].Len(),2) + pow(nbrX[1].Len(),2) + (nbrX[0] * nbrX[1])) * (coorT * cross_NbrX) / 60;
		}
	}

	return 0;
}

int Smooth_BDYODT_FeaturePreserving(const int verIndex)
{
	verElement* curVer = verList[verIndex];
	
	double A,B;
	DBL3DVECT dirMainCurvature;
	
	// to compute the direction along which the vertex should move
	double eigenValues[3];
	EIGENVECT possibleDirs = GetEigenVector(curVer,eigenValues);
	DBL3DVECT dir1,dir2;
	dir1.x = possibleDirs.x1;
	dir1.y = possibleDirs.y1;
	dir1.z = possibleDirs.z1;
	dir2.x = possibleDirs.x2;
	dir2.y = possibleDirs.y2;
	dir2.z = possibleDirs.z2;
	DBL3DVECT dirMainCurv;
	dirMainCurv = dir1 / dir2;
	dirMainCurv /= dirMainCurv.Len();
	
	if (0 != getBDY_ODT_FeaturePres_Paras(verIndex,A,B,dirMainCurvature))
	{
		return 0;
	}
	
	double m = -B / (2 * A);
	
	
	curVer->m_x += m * dirMainCurvature.x;
	curVer->m_y += m * dirMainCurvature.y;
	curVer->m_z += m * dirMainCurvature.z;
	
	return 0;
}

int getBDY_ODT_FeaturePres_Paras(const int verIndex,double& A,double& B,const DBL3DVECT& dirMainCurv)
{
	verElement* curVer = verList[verIndex];

	DBL3DVECT curVDBL;
	curVDBL = VerElement2DBLVECT(curVer);

	double triArea;
	DBL3DVECT triNormal;
	double norMag;
	DBL3DVECT normalVector;

	int numTet = curVer->m_numTet;
	int* adjTet = curVer->m_tets;
	int idxCurTet;
	tetElement* curTet;
	int index_NbrsIn_Tet[3];
	DBL3DVECT nbrVertex[3];// the three neighbors of curV in every tet
	DBL3DVECT nbrDifVector[3];// the three difference vectors from curVer to its three tet neighbors
	double volOmega_star;

	triElement* curTri;
	int index_NbrsIn_Tri[2];
	DBL3DVECT	nbrX[2];// the two DIFFERENCE vectors of from curV to its boundary adjacent neighbors. conunter-clockwise 
	DBL3DVECT   oppoNbr;// the opposite vertex of current vertex in a tet which has a face consisted of current vertex and its two neighbors
	DBL3DVECT   sum_NbrX;
	DBL3DVECT   cross_NbrX;


	// first, compute the volume of Omega_0 (which is also Omega_star)
	// and get the normal vector which determains the tangent plane
	volOmega_star = 0;
	normalVector.x = normalVector.y = normalVector.z = 0;
	int itet = 0;
	for (itet = 0; itet < numTet; itet ++)
	{
		curTet = &(tetMESH[curVer->m_tets[itet]]);
		index_NbrsIn_Tet[0] = curVer->m_subscriptIn_tet_2nd[itet];
		index_NbrsIn_Tet[1] = curVer->m_subscriptIn_tet_3rd[itet];
		index_NbrsIn_Tet[2] = curVer->m_subscriptIn_tet_4th[itet];		
		index_NbrsIn_Tet[0] = curTet->m_v[index_NbrsIn_Tet[0]];
		index_NbrsIn_Tet[1] = curTet->m_v[index_NbrsIn_Tet[1]];
		index_NbrsIn_Tet[2] = curTet->m_v[index_NbrsIn_Tet[2]];		

		nbrVertex[0] = VerElement2DBLVECT(verList[index_NbrsIn_Tet[0]]);
		nbrVertex[1] = VerElement2DBLVECT(verList[index_NbrsIn_Tet[1]]);
		nbrVertex[2] = VerElement2DBLVECT(verList[index_NbrsIn_Tet[2]]);		
		volOmega_star += Volume(curVDBL,nbrVertex[0],nbrVertex[1],nbrVertex[2]);


		nbrDifVector[0] = nbrVertex[0] - curVDBL;		

		triArea = Area(nbrVertex[0],nbrVertex[1],nbrVertex[2]);
		triNormal = (nbrVertex[1] - nbrVertex[0]) / (nbrVertex[2] - nbrVertex[0]);
		norMag = triNormal.Len();
		triNormal /= norMag;

		if (triNormal * nbrDifVector[0] > 1e-10)
		{
			triNormal = -triNormal;			
		}
		normalVector = normalVector + triNormal * triArea;
	}

	norMag = normalVector.Len();
	normalVector /= norMag;

	// third, compute the A
	A = volOmega_star / 4;	
 	for (itet = 0; itet < numTet; itet ++)
	{
		idxCurTet = curVer->m_tets[itet];
		curTet = &(tetMESH[idxCurTet]);
		for (int iFace = 0; iFace < 4; iFace ++)
		{
			curTri = &(tetFACES[curTet->m_faces[iFace]]);
			if (curTri->m_type != FACETYPE_BDY)
				continue;
			if (verIndex != curTri->m_v[0] && verIndex != curTri->m_v[1] && verIndex != curTri->m_v[2])
				continue;

			if (verIndex == curTri->m_v[0])
			{
				index_NbrsIn_Tri[0] = curTri->m_v[1];
				index_NbrsIn_Tri[1] = curTri->m_v[2];
			}
			else if (verIndex == curTri->m_v[1])
			{
				index_NbrsIn_Tri[0] = curTri->m_v[0];
				index_NbrsIn_Tri[1] = curTri->m_v[2];
			}
			else if (verIndex == curTri->m_v[2])
			{
				index_NbrsIn_Tri[0] = curTri->m_v[0];
				index_NbrsIn_Tri[1] = curTri->m_v[1];
			}
			else 
			{
				perror("Error happened in getBDY_ODT_FeaturePres_Paras");
			}
			nbrX[0] = VerElement2DBLVECT(verList[index_NbrsIn_Tri[0]]);
			nbrX[1] = VerElement2DBLVECT(verList[index_NbrsIn_Tri[1]]);
			nbrX[0] = nbrX[0] - curVDBL;
			nbrX[1] = nbrX[1] - curVDBL;

			sum_NbrX = nbrX[0] + nbrX[1];		
			cross_NbrX = nbrX[0] / nbrX[1];

			if (idxCurTet == curTri->m_tet[0])
			{
				oppoNbr = VerElement2DBLVECT(verList[curTet->m_v[curTri->m_ssb[0]]]);
			}
			else
			{
				perror("Error");
			}
			oppoNbr = oppoNbr - curVDBL;
		
			if (cross_NbrX * oppoNbr > 1e-10)
			{
				cross_NbrX = -cross_NbrX;
			}


			A -= (dirMainCurv * sum_NbrX) * (dirMainCurv * cross_NbrX) / 60;		
		}
	}	

	if (A > 0)
	{		
	}
	else
	{
 		perror("A < 0!");		

		return -1;
	}	

	// fourth, compute B 
	B = 0;
	for (itet = 0; itet < numTet; itet ++)
	{
		curTet = &(tetMESH[curVer->m_tets[itet]]);
		index_NbrsIn_Tet[0] = curVer->m_subscriptIn_tet_2nd[itet];
		index_NbrsIn_Tet[1] = curVer->m_subscriptIn_tet_3rd[itet];
		index_NbrsIn_Tet[2] = curVer->m_subscriptIn_tet_4th[itet];		
		index_NbrsIn_Tet[0] = curTet->m_v[index_NbrsIn_Tet[0]];
		index_NbrsIn_Tet[1] = curTet->m_v[index_NbrsIn_Tet[1]];
		index_NbrsIn_Tet[2] = curTet->m_v[index_NbrsIn_Tet[2]];		
		
		nbrVertex[0] = VerElement2DBLVECT(verList[index_NbrsIn_Tet[0]]);
		nbrVertex[1] = VerElement2DBLVECT(verList[index_NbrsIn_Tet[1]]);
		nbrVertex[2] = VerElement2DBLVECT(verList[index_NbrsIn_Tet[2]]);		


		nbrDifVector[0] = nbrVertex[0] - curVDBL;		

		triArea = Area(nbrVertex[0],nbrVertex[1],nbrVertex[2]);
		triNormal = (nbrVertex[1] - nbrVertex[0]) / (nbrVertex[2] - nbrVertex[0]);
		norMag = triNormal.Len();
		triNormal /= norMag;
		
		if (triNormal * nbrDifVector[0] > 1e-10)
		{
			triNormal = -triNormal;			
		}

		B += triArea * (dirMainCurv * triNormal) * (pow((nbrVertex[0] - curVDBL).Len(),2) + \
													  pow((nbrVertex[1] - curVDBL).Len(),2) + 
													  pow((nbrVertex[2] - curVDBL).Len(),2));		
	}

	B /= 12;	

	for (itet = 0; itet < numTet; itet ++)
	{
		idxCurTet = curVer->m_tets[itet];
		curTet = &(tetMESH[idxCurTet]);
		for (int iFace = 0; iFace < 4; iFace ++)
		{
			curTri = &(tetFACES[curTet->m_faces[iFace]]);
			if (curTri->m_type != FACETYPE_BDY)
				continue;
			if (verIndex != curTri->m_v[0] && verIndex != curTri->m_v[1] && verIndex != curTri->m_v[2])
				continue;

			if (verIndex == curTri->m_v[0])
			{
				index_NbrsIn_Tri[0] = curTri->m_v[1];
				index_NbrsIn_Tri[1] = curTri->m_v[2];
			}
			else if (verIndex == curTri->m_v[1])
			{
				index_NbrsIn_Tri[0] = curTri->m_v[0];
				index_NbrsIn_Tri[1] = curTri->m_v[2];
			}
			else if (verIndex == curTri->m_v[2])
			{
				index_NbrsIn_Tri[0] = curTri->m_v[0];
				index_NbrsIn_Tri[1] = curTri->m_v[1];
			}
			else
			{
				perror("Error.");
			}
		
			nbrX[0] = VerElement2DBLVECT(verList[index_NbrsIn_Tri[0]]);
			nbrX[1] = VerElement2DBLVECT(verList[index_NbrsIn_Tri[1]]);
			nbrX[0] = nbrX[0] - curVDBL;
			nbrX[1] = nbrX[1] - curVDBL;
			cross_NbrX = nbrX[0] / nbrX[1];
		
			if (idxCurTet == curTri->m_tet[0])
			{
				oppoNbr = VerElement2DBLVECT(verList[curTet->m_v[curTri->m_ssb[0]]]);
			}
			else
			{
				perror("Error.");
			}
			oppoNbr = oppoNbr - curVDBL;
			if (cross_NbrX * oppoNbr  > 1e-10)
			{
				cross_NbrX = -cross_NbrX;				
			}

			B -= (pow(nbrX[0].Len(),2) + pow(nbrX[1].Len(),2) + (nbrX[0] * nbrX[1])) * (dirMainCurv * cross_NbrX) / 60;		
		}
	}

	return 0;
}

EIGENVECT GetEigenVector(verElement* curVer,double eigenValues[3])
{	
	double A[3][3];
	
	for (int iRow = 0; iRow < 3; iRow ++)
	{
		for (int iCol = 0; iCol < 3; iCol ++)
		{
			A[iRow][iCol] = 0;
		}
	}
	
	tetElement* curTet;
	triElement* curFace;
	verElement* tetVers[4];
	double tmpFaceNormal[3];
	double faceArea;
	for (int iTet = 0; iTet < curVer->m_numTet; iTet ++)
	{
		curTet = &(tetMESH[curVer->m_tets[iTet]]);
		for (int iFace = 0; iFace < 4; iFace ++)
		{
			curFace = &(tetFACES[curTet->m_faces[iFace]]);
			if (curFace->m_type != FACETYPE_BDY)
				continue;
			
			tetVers[0] = verList[curFace->m_v[0]];
			tetVers[1] = verList[curFace->m_v[1]];
			tetVers[2] = verList[curFace->m_v[2]];
			tetVers[3] = verList[curTet->m_v[iFace]];			
			ComputeFaceNormal(tetVers[0],tetVers[1],tetVers[2],tetVers[3],tmpFaceNormal,false);
			faceArea = Area(tetVers[0],tetVers[1],tetVers[2]);
			
			for (int iRow = 0; iRow < 3; iRow ++)
			{
				for (int iCol = 0; iCol < 3; iCol ++)
				{
					A[iRow][iCol] += tmpFaceNormal[iRow] * tmpFaceNormal[iCol] * faceArea;
				}
			}
		}			
	}
	// end of construction of A
	
  return EigenVectors_Jocobi(A,eigenValues);  
}

void ComputeFaceNormal( const verElement* const Ver0,\
						  const verElement* const Ver1,\
						  const verElement* const Ver2,\
						  const verElement* const Ver3,\
						  double faceNormal[3],bool PointInside)
{
	DBL3DVECT Dir[2];
	Dir[0].x = Ver1->m_x - Ver0->m_x;
	Dir[0].y = Ver1->m_y - Ver0->m_y;
	Dir[0].z = Ver1->m_z - Ver0->m_z;
	Dir[1].x = Ver2->m_x - Ver0->m_x;
	Dir[1].y = Ver2->m_y - Ver0->m_y;
	Dir[1].z = Ver2->m_z - Ver0->m_z;
	
	DBL3DVECT fNormal = Dir[0] / Dir[1];
	double norMag = sqrt(fNormal.x * fNormal.x + fNormal.y * fNormal.y + fNormal.z * fNormal.z);
	if (fabs(norMag) < 1e-10)
	{
		perror("There are degenerated triangles.");
		return;
	}
	fNormal /= norMag;

	
	DBL3DVECT refDir;
	refDir.x = Ver3->m_x - Ver0->m_x;
	refDir.y = Ver3->m_y - Ver0->m_y;
	refDir.z = Ver3->m_z - Ver0->m_z;
	refDir = refDir / refDir.Len();
	
	//
	faceNormal[0] = fNormal.x;
	faceNormal[1] = fNormal.y;
	faceNormal[2] = fNormal.z;
	
	if (fNormal * refDir < 0)
	{
		if (PointInside)
		{
			faceNormal[0] = -faceNormal[0];
			faceNormal[1] = -faceNormal[1];
			faceNormal[2] = -faceNormal[2];
		}
	}
	else
	{
		if (!PointInside)
		{
			faceNormal[0] = -faceNormal[0];
			faceNormal[1] = -faceNormal[1];
			faceNormal[2] = -faceNormal[2];
		}
	}
}

EIGENVECT EigenVectors_Jocobi(double Matrix[3][3],double eigenValues[3])
{
	double X[3];
	double Y[3];
	double Z[3];

	double P[3][3] = {{1,0,0},{0,1,0},{0,0,1}},Q[3][3],classT[3][3];
	double a_rp,a_rq,a_pp,a_qq,a_pq;
	double theta,t,gama;
	double C,S;
	double tao = 0;
	int r = 0, p = 0, q = 0;
	double maxPQEle;
	int row = 0,col = 0;
	int k = 0;	
	for (row = 0; row < 3; row ++)
	{
		for (col = row + 1; col < 3; col ++)
		{				
			tao += Matrix[row][col] * Matrix[row][col];
		}
	}
	
	while(tao > 1e-10)
	{
		// find the max pq element
		maxPQEle = fabs(Matrix[0][1]);
		p = 0; q = 1;
		for(row = 0; row < 3; row ++)
		{
			for (col = row + 1; col < 3; col ++)
			{
				if (fabs(Matrix[row][col]) > maxPQEle)
				{
					maxPQEle = fabs(Matrix[row][col]);
					p = row; q = col;
				}
			}
		}
		r = 3 - p - q;
		
		// eliminate the p,q elements
		if (fabs(Matrix[p][q]) < 1e-10)
			break;
		else
		{
			theta = (Matrix[q][q] - Matrix[p][p]) / (2 * Matrix[p][q]);
			if (1.0 / fabs(theta) < 1e-10)
				t = 1.0 / (2 * theta);
			else
				t = theta / fabs(theta) / (fabs(theta) + sqrt(theta * theta + 1));
			C = 1.0 / sqrt(t * t + 1);
			S = t * C;
			gama = S / (1 + C);
			a_rp = Matrix[r][p]; a_rq = Matrix[r][q];
			a_pp = Matrix[p][p]; a_qq = Matrix[q][q];
			a_pq = Matrix[p][q];
			Matrix[p][p] = a_pp - t * a_pq;
			Matrix[q][q] = a_qq + t * a_pq;
			Matrix[r][p] = Matrix[p][r] = a_rp - S * (a_rq + gama * a_rp);
			Matrix[r][q] = Matrix[q][r] = a_rq + S * (a_rp - gama * a_rq);
			Matrix[p][q] = Matrix[q][p] = (C * C - S * S) * a_pq + C * S * (a_pp - a_qq);
			
			// compute the P matrix
			for (row = 0; row < 3; row ++)
			{
				for (col = 0; col < 3; col ++)
				{
					Q[row][col] = P[row][col];
					classT[row][col] = 0;
				}
				classT[row][row] = 1;
			}
			classT[p][p] = C; classT[q][q] = C;
			classT[p][q] = S; classT[q][p] = -S;
			for (row = 0; row < 3; row ++)
			{
				for (col = 0; col < 3; col ++)
				{
					P[row][col] = 0;
					for (k = 0; k < 3; k ++)
					{
						P[row][col] += Q[row][k] * classT[k][col];
					}
				}
			}
			
		}
		
		// recompute the tao value.
		tao = 0;
		for (row = 0; row < 3; row ++)
		{
			for (col = row + 1; col < 3; col ++)
			{				
				tao += Matrix[row][col] * Matrix[row][col];
			}
		}
	}
	
	// sort the eigenvalues and reordering the eigenvectors. Suppose the index of eigenvalues in decreasing order is p,q,r 
	// find the max eigenvalue and interchange the elements of Matrix and P
	maxPQEle = Matrix[0][0];p = 0;
	for (row = 1; row < 3; row ++)
	{
		if (Matrix[row][row] > maxPQEle)
		{
			maxPQEle = Matrix[row][row];
			p = row;
		}
	}
	double temp = Matrix[0][0];
	Matrix[0][0] = Matrix[p][p];
	Matrix[p][p] = temp;
	for (row = 0; row < 3; row ++)
	{
		temp = P[row][0];
		P[row][0] = P[row][p];
		P[row][p] = temp;
	}
	
	// find the min eigenvalue and interchange the elements of Matrix and P
	maxPQEle = Matrix[1][1];
	if (maxPQEle < Matrix[2][2])
	{// interchange the 1 and 2 elements
		temp = Matrix[1][1];
		Matrix[1][1] = Matrix[2][2];
		Matrix[2][2] = temp;
		
		for (row = 0; row < 3; row ++)
		{
			temp = P[row][1];
			P[row][1] = P[row][2];
			P[row][2] = temp;
		}
	}
	
	for (row = 0; row < 3; row ++)
	{
		Z[row] = P[row][0];// corresponding to the max eigenvalue
		X[row] = P[row][1];// corresponding to the medium eigenvalue
		Y[row] = P[row][2];// corresponding to the min eigenvalue
	}
	temp = sqrt(Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2]);
	Z[0] /= temp; Z[1] /= temp; Z[2] /= temp;
	temp = sqrt(Y[0] * Y[0] + Y[1] * Y[1] + Y[2] * Y[2]);
	Y[0] /= temp; Y[1] /= temp; Y[2] /= temp;
	temp = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
	X[0] /= temp; X[1] /= temp; X[2] /= temp;
	
	EIGENVECT eigenVectors;
	eigenVectors.x1 = Z[0];
	eigenVectors.y1 = Z[1];
	eigenVectors.z1 = Z[2];
	eigenVectors.x2 = X[0];
	eigenVectors.y2 = X[1];
	eigenVectors.z2 = X[2];
	eigenVectors.x3 = Y[0];
	eigenVectors.y3 = Y[1];
	eigenVectors.z3 = Y[2];

	eigenValues[0] = Matrix[0][0];
	eigenValues[1] = Matrix[1][1];
	eigenValues[2] = Matrix[2][2];

	return eigenVectors;
}

DBL3DVECT VerElement2DBLVECT(const verElement* ver)
{
	DBL3DVECT tmp;
	tmp.x = ver->m_x;
	tmp.y = ver->m_y;
	tmp.z = ver->m_z;
	
	return tmp;
}

void DBLVECT2VerElement(const DBL3DVECT &DBL, verElement *&ver)
{
	ver->m_x = DBL.x;
	ver->m_y = DBL.y;
	ver->m_z = DBL.z;
}

int ComputeDihedralAngles_All(void)
{	
	if (tetMESH.size() == 0)
		return -1;

	double localMinAngle = 0,localMaxAngle = 0;
	double minAngle = 181,maxAngle = -1;
	
	for (int itet = 0; itet < m_numTet; itet ++)
	{
		ComputeDiAngle_SigleTet(&(tetMESH[itet]), localMinAngle, localMaxAngle);
		
		minAngle = min(minAngle,localMinAngle);
		maxAngle = max(maxAngle,localMaxAngle);
	}

	m_minMinAngle = minAngle;
	m_maxMaxAngle = maxAngle;
	
	return 0;
}

double ComputeLocalAngle(const int verIndex,double& localMaxAngle)
{
	double minDiAngle = 181;
	double maxDiAngle = -1;
	
	const verElement* curVer = verList[verIndex];
	const int numAdjTets = curVer->m_numTet;
	const int* adjTets = curVer->m_tets;	
	
	const DBL3DVECT curVerCoord = VerElement2DBLVECT(curVer);

	int kTet;
	const tetElement* curTet = NULL;
	DBL3DVECT faceNodes[3];
	double curAngle = 0;
	for (kTet = 0; kTet < numAdjTets; kTet ++)
	{
		curTet = &(tetMESH[adjTets[kTet]]);
		
		faceNodes[0] = VerElement2DBLVECT(verList[curTet->m_v[(curVer->m_subscriptIn_tet[kTet] + 1) % 4]]);
		faceNodes[1] = VerElement2DBLVECT(verList[curTet->m_v[(curVer->m_subscriptIn_tet[kTet] + 2) % 4]]);
		faceNodes[2] = VerElement2DBLVECT(verList[curTet->m_v[(curVer->m_subscriptIn_tet[kTet] + 3) % 4]]);

		curAngle = ComputeMinDiAngle(curVerCoord,faceNodes[0],faceNodes[1],faceNodes[2]);
		
		minDiAngle = min(minDiAngle,curAngle);		
		maxDiAngle = max(maxDiAngle,curAngle);
	}
	
	return minDiAngle;
}

double ComputeMinDiAngle(const DBL3DVECT& P1,const DBL3DVECT& P2,const DBL3DVECT& P3,const DBL3DVECT& P4)
{
	double minAngle;
	
	double curDiAngle;
	minAngle = ComputeDiAngle(P1,P2,P3,P4);	
	
	curDiAngle = ComputeDiAngle(P1,P3, P2,P4);
	minAngle = min(minAngle,curDiAngle); 	
	
	curDiAngle = ComputeDiAngle(P1,P4, P2,P3);
	minAngle = min(minAngle,curDiAngle); 
	
	curDiAngle = ComputeDiAngle(P2,P3, P1,P4);
	minAngle = min(minAngle,curDiAngle); 
	
	curDiAngle = ComputeDiAngle(P2,P4, P1,P3);
	minAngle = min(minAngle,curDiAngle); 
	
	curDiAngle = ComputeDiAngle(P3,P4, P1,P2);
	minAngle = min(minAngle,curDiAngle); 
	
	return minAngle;
}

int ComputeDiAngle_SigleTet(const tetElement* curTet, double& minAngle,double& maxAngle)
{
	verElement* verPt[4] = {verList[curTet->m_v[0]],\
							verList[curTet->m_v[1]],\
							verList[curTet->m_v[2]],\
							verList[curTet->m_v[3]]};
	double curDiAngle;
	minAngle = maxAngle = ComputeDiAngle(verPt[0],verPt[1],verPt[2],verPt[3]);
	
	curDiAngle = ComputeDiAngle(verPt[0],verPt[2], verPt[1],verPt[3]);
	minAngle = min(minAngle,curDiAngle);
	maxAngle = max(maxAngle, curDiAngle);
	
	curDiAngle = ComputeDiAngle(verPt[0],verPt[3], verPt[1],verPt[2]);
	minAngle = min(minAngle,curDiAngle); 
	maxAngle = max(maxAngle, curDiAngle);
	
	curDiAngle = ComputeDiAngle(verPt[1],verPt[2], verPt[0],verPt[3]);
	minAngle = min(minAngle,curDiAngle); 
	maxAngle = max(maxAngle, curDiAngle);
	
	curDiAngle = ComputeDiAngle(verPt[1],verPt[3], verPt[0],verPt[2]);
	minAngle = min(minAngle,curDiAngle); 
	maxAngle = max(maxAngle, curDiAngle);
	
	curDiAngle = ComputeDiAngle(verPt[2],verPt[3], verPt[0],verPt[1]);
	minAngle = min(minAngle,curDiAngle); 
	maxAngle = max(maxAngle, curDiAngle);
	
	return 0;
}

double ComputeDiAngle(const verElement* P,const verElement* Q,const verElement* r1,const verElement* r2)
{// PQ is the axis to compute dihedral angle
	DBL3DVECT tmpP;
	DBL3DVECT tmpQ;
	DBL3DVECT tmpR1;
	DBL3DVECT tmpR2;
	tmpP = VerElement2DBLVECT(P);
	tmpQ = VerElement2DBLVECT(Q);
	tmpR1 = VerElement2DBLVECT(r1);
	tmpR2 = VerElement2DBLVECT(r2);
	return ComputeDiAngle(tmpP,tmpQ,tmpR1,tmpR2);
}

double ComputeDiAngle(const DBL3DVECT& P,const DBL3DVECT& Q,const DBL3DVECT& r1,const DBL3DVECT& r2)
{
	DBL3DVECT normPQ1 = (Q - P) / (r1 - P);
	DBL3DVECT normPQ2 = (Q - P) / (r2 - P);// normP12 and normQ21 has the same perpoty. e.g. 
									// they point to the inner or outer side of the tethedral simultaneously.
	DBL3DVECT Pr2;
	Pr2.x = r2.x - P.x;
	Pr2.y = r2.y - P.y;
	Pr2.z = r2.z - P.z;
	DBL3DVECT Pr1;
	Pr1.x = r1.x - P.x;
	Pr1.y = r1.y - P.y;
	Pr1.z = r1.z - P.z;
	if ((normPQ1 * Pr2) < 0)
	{
		normPQ1.x = -normPQ1.x;
		normPQ1.y = -normPQ1.y;
		normPQ1.z = -normPQ1.z;
	}
	if ((normPQ2 * Pr1) < 0)
	{
		normPQ2.x = -normPQ2.x;
		normPQ2.y = -normPQ2.y;
		normPQ2.z = -normPQ2.z;
	}
	
	double dotPro = normPQ1 * normPQ2;
	double nor1 = normPQ1.Len();
	double nor2 = normPQ2.Len();
	double cosValue = dotPro / nor1 / nor2;
	double angle;
	if (fabs(cosValue - 1) < 1e-10)
		angle = 0;
	else if (fabs(cosValue + 1) < 1e-10)
		angle = MATHPI;
	else
		angle = acos(cosValue);
	
	return 180 - angle / MATHPI * 180;
}

bool CheckValid_byPosVol(const int verIndex)
{
	const verElement* curVer = verList[verIndex];
	
	const int numTet = curVer->m_numTet;
	const tetElement* curTet = NULL;
	double tetSign;
	for (int kTet = 0; kTet < numTet; kTet ++)
	{
		curTet = &(tetMESH[curVer->m_tets[kTet]]);
		tetSign = SignedVolume(curTet);
		
		if (tetSign * m_VolumeSignes[curVer->m_tets[kTet]] < 0) 
		{
			return false;
		}
	}
	
	return true;
}

double Area(const DBL3DVECT& P1,const DBL3DVECT& P2,const DBL3DVECT& P3)
{// the area of the triangle formed by P1,P2,P3
	return ((P2 - P1) / (P3 - P1)).Len() / 2;
}

double Area(const verElement* v1, const verElement* v2,const verElement* v3)
{
	DBL3DVECT P1,P2,P3;
	
	P1 = VerElement2DBLVECT(v1);
	P2 = VerElement2DBLVECT(v2);
	P3 = VerElement2DBLVECT(v3);
	
	return Area(P1,P2,P3);
}

double Volume(const DBL3DVECT& P1,const DBL3DVECT& P2,const DBL3DVECT& P3,const DBL3DVECT& P4)
{// volume of the tet formed by P1 to P4
	DBL3DVECT DP2,DP3,DP4;
	DP2.x = P2.x - P1.x; DP2.y = P2.y - P1.y; DP2.z = P2.z - P1.z;
	DP3.x = P3.x - P1.x; DP3.y = P3.y - P1.y; DP3.z = P3.z - P1.z;
	DP4.x = P4.x - P1.x; DP4.y = P4.y - P1.y; DP4.z = P4.z - P1.z;
	
	DBL3DVECT theCross = DP2 / DP3;	
	
	double theVolume = fabs(DP4 * theCross);
	
	return theVolume / 6.0;
}

double SignedVolume(const DBL3DVECT& P1,const DBL3DVECT& P2,const DBL3DVECT& P3,const DBL3DVECT& P4)
{
	DBL3DVECT DP2,DP3,DP4;
	DP2.x = P2.x - P1.x; DP2.y = P2.y - P1.y; DP2.z = P2.z - P1.z;
	DP3.x = P3.x - P1.x; DP3.y = P3.y - P1.y; DP3.z = P3.z - P1.z;
	DP4.x = P4.x - P1.x; DP4.y = P4.y - P1.y; DP4.z = P4.z - P1.z;
	
	DBL3DVECT theCross = DP2 / DP3;	
	
	double theVolume = DP4 * theCross;
	
	return theVolume / 6.0;
}

double SignedVolume(const tetElement* curTet)
{
	DBL3DVECT Vers[4];
	Vers[0] = VerElement2DBLVECT(verList[curTet->m_v[0]]);
	Vers[1] = VerElement2DBLVECT(verList[curTet->m_v[1]]);
	Vers[2] = VerElement2DBLVECT(verList[curTet->m_v[2]]);
	Vers[3] = VerElement2DBLVECT(verList[curTet->m_v[3]]);
	
	return SignedVolume(Vers[0],Vers[1],Vers[2],Vers[3]);
}

void InitVolumeSigns()
{
	const tetElement* curTet = NULL;;
	
	for (int kTet = 0; kTet < m_numTet; kTet ++)
	{
		curTet = &(tetMESH[kTet]);
		
		(SignedVolume(curTet) > 0) ? (m_VolumeSignes[kTet] = 1) : (m_VolumeSignes[kTet] = -1);	
	}
}

int InitBasicData(void)
{	
	m_numVer = m_numTet = m_numTetFACES /*= m_numTri*/ = 0;
	
	ori_minX = ori_maxX = ori_minY = ori_maxY = ori_minZ = ori_maxZ = 0;	
	
	m_oriMetric = 0;

	return 0;
}	

int LoadData(const char* filename)
{
	int meshType = 0;
	double x = 0,y = 0,z = 0;	
	FILE* fp = NULL;	

	fp = fopen(filename,"r");
	if (fp == NULL)
	{
		perror("Open Input File Error.");
		return -1;
	}
	ReleaseAll();
	InitBasicData();
	
	char fileFormat[32];
	fscanf(fp,"%s",fileFormat);
	if (strcmp(fileFormat,"TET3") != 0)
	{
		perror("File format is not TET3.");
		return -1;
	}
	
	triElement newFace;
	fscanf(fp,"%d%d%d",&m_numVer,&m_numTet,&m_numTetFACES);
	tetFACES.assign(m_numTetFACES,newFace);
	
	tetElement newTet;
	try{
		verList.assign(m_numVer,NULL);
		tetMESH.assign(m_numTet,newTet);		
	}
	catch(...)
	{
		perror("Not enough memory.");
		return -1;
	}

	verElement* newVer = new verElement[1];
	fscanf(fp,"%lf%lf%lf",&x,&y,&z);
	ori_minX = ori_maxX = newVer->m_x = x;
	ori_minY = ori_maxY = newVer->m_y = y;
	ori_minZ = ori_maxZ = newVer->m_z = z;
	
	verList[0] = newVer;
	int iver = 0;
	for (iver = 1; iver < m_numVer; iver ++)
	{
 		newVer = new verElement[1];
		fscanf(fp,"%lf%lf%lf",&x,&y,&z);
		newVer->m_x = x; newVer->m_y = y; newVer->m_z = z;
		verList[iver] = newVer;

			ori_minX = min(ori_minX,x); ori_maxX = max(ori_maxX,x); 
			ori_minY = min(ori_minY,y); ori_maxY = max(ori_maxY,y);
			ori_minZ = min(ori_minZ,z); ori_maxZ = max(ori_maxZ,z);
	}
	
	m_oriMetric = (ori_maxX - ori_minX + ori_maxY - ori_minY + ori_maxZ - ori_minZ) / 3;

	triElement* tmpFace;
	tetElement* tmpTet;
	int itet = 0,iface = 0;
	for (itet = 0; itet < m_numTet; itet ++)
	{
		fscanf(fp,"%d %d%d%d%d %d%d%d%d  %d%d %d%d %d%d %d%d",&meshType,\
					&(tetMESH[itet].m_v[0]),&(tetMESH[itet].m_v[1]),&(tetMESH[itet].m_v[2]),&(tetMESH[itet].m_v[3]),\
					&(tetMESH[itet].m_triType[0]),&(tetMESH[itet].m_triType[1]),&(tetMESH[itet].m_triType[2]),&(tetMESH[itet].m_triType[3]),\
					&(tetMESH[itet].m_tets[0]),&(tetMESH[itet].m_faces[0]),\
					&(tetMESH[itet].m_tets[1]),&(tetMESH[itet].m_faces[1]),\
					&(tetMESH[itet].m_tets[2]),&(tetMESH[itet].m_faces[2]),\
					&(tetMESH[itet].m_tets[3]),&(tetMESH[itet].m_faces[3]));

				tmpTet = &(tetMESH[itet]);				
	}
	for (iface = 0; iface < m_numTetFACES; iface ++)
	{
		fscanf(fp,"%d%d%d %d%d %d%d",&(tetFACES[iface].m_v[0]),&(tetFACES[iface].m_v[1]),&(tetFACES[iface].m_v[2]),\
									&(tetFACES[iface].m_tet[0]),&(tetFACES[iface].m_ssb[0]),\
									&(tetFACES[iface].m_tet[1]),&(tetFACES[iface].m_ssb[1]));

		tmpFace = &(tetFACES[iface]);				
	}
	fclose(fp);

	ComputeTetInfo_forVer();
		
	ComputeBdyInfo_forTet();

	CorrectTetInfo();
	ComputeAllFaceNormal();
	
	ComputeBdyVerNormal();

	ComputeDihedralAngles_All();

	// statistic the number of BOUNDARY faces
	int faceNum = 0,kFace = 0;
	const int numTetFaces = tetFACES.size();
	for (kFace = 0; kFace < numTetFaces; kFace ++)
	{
		if (tetFACES[kFace].m_tet[1] == -1)
			faceNum ++;
	}
		
	try{
		m_VolumeSignes.assign(m_numTet,0);		
	}
	catch(...)
	{
		perror("Not enough memory for 'm_InitVolumeSignes'.\n");
	}
	InitVolumeSigns();

	return 0;
}

int SaveData(const char* filename)
{
	if (verList.size() == 0 || tetMESH.size() == 0)
		return -1;

	FILE* fp = NULL;
	fp = fopen(filename,"w");
	if (fp == NULL)
	{
		perror("File open error!");
		return -1;
	}
	fprintf(fp,"TET3\n");
	fprintf(fp,"%d %d %d\n",m_numVer,m_numTet,m_numTetFACES);
	
	int iver = 0;
	for (iver = 0; iver < m_numVer; iver ++)
	{
		fprintf(fp,"%lf %lf %lf\n",verList[iver]->m_x,verList[iver]->m_y,verList[iver]->m_z);
	}

	int itet = 0;
	for (itet = 0; itet < m_numTet; itet ++)
	{
		fprintf(fp,"4 %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",\
					tetMESH[itet].m_v[0],tetMESH[itet].m_v[1],tetMESH[itet].m_v[2],tetMESH[itet].m_v[3],\
					tetMESH[itet].m_triType[0],tetMESH[itet].m_triType[1],tetMESH[itet].m_triType[2],tetMESH[itet].m_triType[3],\
					tetMESH[itet].m_tets[0],tetMESH[itet].m_faces[0],\
					tetMESH[itet].m_tets[1],tetMESH[itet].m_faces[1],\
					tetMESH[itet].m_tets[2],tetMESH[itet].m_faces[2],\
					tetMESH[itet].m_tets[3],tetMESH[itet].m_faces[3]);
	}

	int iface = 0;
	for (iface = 0; iface < m_numTetFACES; iface ++)
	{
		fprintf(fp,"%d %d %d %d %d %d %d\n",\
					tetFACES[iface].m_v[0],tetFACES[iface].m_v[1],tetFACES[iface].m_v[2],\
					tetFACES[iface].m_tet[0],tetFACES[iface].m_ssb[0],\
					tetFACES[iface].m_tet[1],tetFACES[iface].m_ssb[1]);
	}
	
	fclose(fp);

	return 0;
}

int ComputeTetInfo_forVer(void)
{
	if (verList.size() == 0 || tetMESH.size() == 0)
		return -1;

	std::vector<int>* tmpTetList = new std::vector<int>[m_numVer];
	std::vector<int>* tmpRelaSubscript = new std::vector<int>[m_numVer];
	std::vector<int>* tmpRelaSubscript2 = new std::vector<int>[m_numVer];
	std::vector<int>* tmpRelaSubscript3 = new std::vector<int>[m_numVer];
	std::vector<int>* tmpRelaSubscript4 = new std::vector<int>[m_numVer];
	
	// traversal all the tet elements
	int* vers;// all the index of vertex in cur tet
	tetElement* curTet;
 	int numTri = 0;
	int itet = 0;
	for (itet = 0; itet < m_numTet; itet ++)
	{// to get the information of adjacent tets for every vertex.
		curTet = &(tetMESH[itet]);
		vers = curTet->m_v;// all the 4 vertex index in current tet
		
		tmpTetList[vers[0]].push_back(itet);// to record that the current tet is adjacent to the vers[0]-th vertex
		tmpTetList[vers[1]].push_back(itet);
		tmpTetList[vers[2]].push_back(itet);
		tmpTetList[vers[3]].push_back(itet);
		
		tmpRelaSubscript[vers[0]].push_back(0);// to record that the vers[0]-th vertex are the 0'th vertex in current tet
		tmpRelaSubscript2[vers[0]].push_back(1);
		tmpRelaSubscript3[vers[0]].push_back(2);
		tmpRelaSubscript4[vers[0]].push_back(3);// to record that the neighbors of the vers[0]-th vertex in current tet is the 1st, 2nd and 3rd vertex in the current tet

		tmpRelaSubscript[vers[1]].push_back(1);
		tmpRelaSubscript2[vers[1]].push_back(0);
		tmpRelaSubscript3[vers[1]].push_back(2);
		tmpRelaSubscript4[vers[1]].push_back(3);

		tmpRelaSubscript[vers[2]].push_back(2);
		tmpRelaSubscript2[vers[2]].push_back(1);
		tmpRelaSubscript3[vers[2]].push_back(0);
		tmpRelaSubscript4[vers[2]].push_back(3);

		tmpRelaSubscript[vers[3]].push_back(3);
		tmpRelaSubscript2[vers[3]].push_back(1);
		tmpRelaSubscript3[vers[3]].push_back(2);
		tmpRelaSubscript4[vers[3]].push_back(0);
		
		if (curTet->m_triType[0] == 0)// 0 is boundary tri
		{
			verList[vers[1]]->m_type = 0;
			verList[vers[2]]->m_type = 0;
			verList[vers[3]]->m_type = 0;// set all the 3 vertex to be boundary vertex			

 			numTri ++;
		}
		if (curTet->m_triType[1] == 0)// 1 is boundary tri
		{
			verList[vers[0]]->m_type = 0;
			verList[vers[2]]->m_type = 0;
			verList[vers[3]]->m_type = 0;
			
 			numTri ++;// Update the number of boundary triangles.
		}
		if (curTet->m_triType[2] == 0)// 2 is boundary tri
		{
			verList[vers[0]]->m_type = 0;
			verList[vers[1]]->m_type = 0;
			verList[vers[3]]->m_type = 0;

 			numTri ++;
		}
		if (curTet->m_triType[3] == 0)// 3 is boundary tri
		{
			verList[vers[0]]->m_type = 0;
			verList[vers[1]]->m_type = 0;
			verList[vers[2]]->m_type = 0;

 			numTri ++;
		}
	}

	// initial all the nodes other than boundary ones as INNER node (23 June 2011)
	for (int kVer = 0; kVer < m_numVer; kVer ++)
	{
		if (verList[kVer]->m_type == -1)
			verList[kVer]->m_type = PTTYPE_INN;
	}
	
	// Update vertex information
	verElement* curVer;// cur vertex
	int tmpNumTet;// the  tet number for every vertex
	std::vector<int>* curTetList;// the tet array for every vertex
	std::vector<int>* curRelatIndex;// the index list of current vertex
	std::vector<int>* RelatIndex2;
	std::vector<int>* RelatIndex3;
	std::vector<int>* RelatIndex4;
	int iver = 0;
	for (iver = 0; iver < m_numVer; iver ++)
	{
		curVer = verList[iver];
		curTetList = tmpTetList + iver;
		curRelatIndex = tmpRelaSubscript + iver;
		RelatIndex2 = tmpRelaSubscript2 + iver;
		RelatIndex3 = tmpRelaSubscript3 + iver;
		RelatIndex4 = tmpRelaSubscript4 + iver;
		
		tmpNumTet = curVer->m_numTet = curTetList->size();// Update the m_numTet for every vertex

		if (tmpNumTet == 0)
		{
			perror("There is isolated vertex in current model");
			curVer->m_type = -1;
			continue;
		}
		
		if (curVer->m_tets != NULL)
			delete[] curVer->m_tets;
		curVer->m_tets = new int[tmpNumTet];
		if (curVer->m_subscriptIn_tet != NULL)
			delete[] curVer->m_subscriptIn_tet;
		curVer->m_subscriptIn_tet = new int[tmpNumTet];
		if (curVer->m_subscriptIn_tet_2nd != NULL)
			delete[] curVer->m_subscriptIn_tet_2nd;
		curVer->m_subscriptIn_tet_2nd = new int[tmpNumTet];
		if (curVer->m_subscriptIn_tet_3rd != NULL)
			delete[] curVer->m_subscriptIn_tet_3rd;
		curVer->m_subscriptIn_tet_3rd = new int[tmpNumTet];
		if (curVer->m_subscriptIn_tet_4th != NULL)
			delete[] curVer->m_subscriptIn_tet_4th;
		curVer->m_subscriptIn_tet_4th = new int[tmpNumTet];
		
		for (itet = 0; itet < tmpNumTet; itet ++)
		{
			curVer->m_tets[itet] = (*curTetList)[itet];// Update the m_tet member
			curVer->m_subscriptIn_tet[itet] = (*curRelatIndex)[itet];// Update the m_subscriptIn_tet member
			curVer->m_subscriptIn_tet_2nd[itet] = (*RelatIndex2)[itet];
			curVer->m_subscriptIn_tet_3rd[itet] = (*RelatIndex3)[itet];
			curVer->m_subscriptIn_tet_4th[itet] = (*RelatIndex4)[itet];
		}
		
		tmpTetList[iver].clear();
		tmpRelaSubscript[iver].clear();
		tmpRelaSubscript2[iver].clear();
		tmpRelaSubscript3[iver].clear();
		tmpRelaSubscript4[iver].clear();
	}
	delete[] tmpTetList;
	delete[] tmpRelaSubscript;
	delete[] tmpRelaSubscript2;
	delete[] tmpRelaSubscript3;
	delete[] tmpRelaSubscript4;

	return 0;
}

int ComputeBdyInfo_forTet(void)
{	
	tetElement* curTet;
	int* vers;// to point the 4 indexes of vertexes of current tet
	for (int itet = 0; itet < m_numTet; itet ++)
	{
		curTet = &(tetMESH[itet]);
		
		// Compute the m_type member
		vers = curTet->m_v;
		if (verList[vers[0]]->m_type * verList[vers[1]]->m_type * verList[vers[2]]->m_type * verList[vers[3]]->m_type == 0)// there are at least 1 vertex on the boundary
			curTet->m_type = 0;
		else
			curTet->m_type = 255;
	}
	
	return 0;
}

int CorrectTetInfo(void)
{
	triElement* curFace;
	for (int iface = 0; iface < m_numTetFACES; iface ++)
	{
		curFace = &(tetFACES[iface]);
		if (curFace->m_tet[0] == -1)
		{
			perror("Error, input face tet index.\n");
		}
		if (curFace->m_tet[1] == -1)
		{
			curFace->m_type = FACETYPE_BDY;
		}
		else
		{
			curFace->m_type = FACETYPE_INN;
		}
	}
	
	return 0;
}

DBL3DVECT ComputeFaceNormal(const triElement* face,bool PointInside)
{
	DBL3DVECT faceVers[3];
	faceVers[0] = VerElement2DBLVECT(verList[face->m_v[0]]);
	faceVers[1] = VerElement2DBLVECT(verList[face->m_v[1]]);
	faceVers[2] = VerElement2DBLVECT(verList[face->m_v[2]]);
	
	DBL3DVECT fNormal = (faceVers[1] - faceVers[0]) / (faceVers[2] - faceVers[0]);
	double norMag = fNormal.Len();
	if (fabs(norMag) < 1e-10)
	{
		perror("There are degenerated triangles.");
		return fNormal;
	}
	fNormal /= norMag;
	
	tetElement* nbTet = &(tetMESH[face->m_tet[0]]);
	DBL3DVECT nbVer;
	nbVer = VerElement2DBLVECT(verList[nbTet->m_v[face->m_ssb[0]]]);
	
	DBL3DVECT refDir;
	refDir = nbVer - faceVers[0]; 
	
	if (fNormal * refDir < 0)
	{
		if (PointInside)
		{
			fNormal = -fNormal;			
		}
	}	
	else
	{
		if (!PointInside)
		{
			fNormal = -fNormal;			
		}
	}
	
	return fNormal;
}

void ComputeAllFaceNormal(void)
{
	triElement* curFace;
	DBL3DVECT faceNormal;
	for (int iFace = 0; iFace < m_numTetFACES; iFace ++)
	{
		curFace = &(tetFACES[iFace]);
		faceNormal = ComputeFaceNormal(curFace);
		curFace->m_nx = faceNormal.x;
		curFace->m_ny = faceNormal.y;
		curFace->m_nz = faceNormal.z;
	}	
}

int ComputeBdyVerNormal(void)
{	
	double curWeight;
	DBL3DVECT faceNormal;
	DBL3DVECT refVect;
	DBL3DVECT verNormal;
	verElement* curVer;
	tetElement* curTet;
	DBL3DVECT NBVers[3];
	for (int iVer = 0; iVer < m_numVer; iVer ++)
	{
		curVer = verList[iVer];
		if (curVer->m_type != PTTYPE_BDY)
			continue;
		
		verNormal = DBL3DVECT(0,0,0);
		for (int iTet = 0; iTet < curVer->m_numTet; iTet ++)
		{
			curTet = &(tetMESH[curVer->m_tets[iTet]]);
			for (int iFace = 0; iFace < 4; iFace ++)
			{
				if (curTet->m_triType[iFace] != FACETYPE_BDY)
					continue;
				if (iVer != curTet->m_v[(iFace + 1) % 4] && iVer != curTet->m_v[(iFace + 2) % 4] && iVer != curTet->m_v[(iFace + 3) % 4])
					continue;
				
				NBVers[0] = VerElement2DBLVECT(verList[curTet->m_v[(iFace + 1) % 4]]);
				NBVers[1] = VerElement2DBLVECT(verList[curTet->m_v[(iFace + 2) % 4]]);
				NBVers[2] = VerElement2DBLVECT(verList[curTet->m_v[(iFace + 3) % 4]]);
				refVect.x = verList[curTet->m_v[iFace]]->m_x - NBVers[0].x;
				refVect.y = verList[curTet->m_v[iFace]]->m_y - NBVers[0].y;
				refVect.z = verList[curTet->m_v[iFace]]->m_z - NBVers[0].z;
				refVect /= refVect.Len();
				
				faceNormal = (NBVers[1] - NBVers[0]) / (NBVers[2] - NBVers[0]);
				faceNormal /= faceNormal.Len();
				
				if (faceNormal * refVect < 0)
				{
					faceNormal = -faceNormal;					
				}
				
				curWeight = Area(NBVers[0],NBVers[1],NBVers[2]);
				verNormal = verNormal + faceNormal * curWeight;
			}
		}
		
		verNormal /= verNormal.Len();
		curVer->m_normal[0] = verNormal.x;
		curVer->m_normal[1] = verNormal.y;
		curVer->m_normal[2] = verNormal.z;
	}
	
	return 0;
}

int ReleaseAll(void)
{	
	int numVer = verList.size();
	for (int iver = 0; iver < numVer; iver ++)
	{
		delete[] verList[iver];
	}
	if (numVer > 0)
		verList.clear();

	if (tetMESH.size() > 0)
		tetMESH.clear();
	if (m_VolumeSignes.size() > 0)
		m_VolumeSignes.clear();

	if (tetFACES.size() > 0)
		tetFACES.clear();

	return 0;
}
