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
#include "tet2tet3.h"

// global variables
int							m_numVer = 0;
std::vector<verElement*>	verList;

int							m_numTetFACES = 0;
std::vector<triElement>		tetFACES;

int							m_numTet = 0;
std::vector<tetElement>		tetMESH;

int Usage(void)
{
	printf("Usage: tet2tet3 <tet2 file>\n");
	return 0;
}

int main(int argc,char** argv)
{
	if (argc < 2)
	{
		Usage();
		return -1;
	}

	const char* inputFileName = argv[1];
	char outputFileName[256] = "\0";
	strcpy(outputFileName,inputFileName);
	const int lenFileName = strlen(outputFileName);
	outputFileName[lenFileName] = '3';
	outputFileName[lenFileName + 1] = '\0';
		
	printf("Initializing basic data....\n");
	InitBasicData();
	
	printf("Loading Data....\n");
	LoadData(inputFileName);

	printf("Reconstructing the boundary information......\n");
	reComputeBdyInfoForTet();

// 	SaveData("mediate");
// 	ReleaseAll();
// 	InitBasicData();
// 	LoadData("mediate");

	printf("Generating the faces of every tetrahedron.....\n");
	GenerateTetFaces();

	printf("Saving result......\n");
	SaveData(outputFileName);


	return 0;
}

int reComputeBdyInfoForTet(void)
{
	if (verList.size() == 0 || tetMESH.size() == 0)
		return -1;
	
	tetElement* curTet;
	int curFaceVerIndex[4];
	int itet = 0;
	for (itet = 0; itet < m_numTet; itet ++)
	{
		((itet % (m_numTet / 10)) == 0) ? (printf("processed %d tet, altogether %d.\n",itet,m_numTet)) : (0);
		curTet = &(tetMESH[itet]);
		curFaceVerIndex[0] = curTet->m_v[0];
		curFaceVerIndex[1] = curTet->m_v[1];
		curFaceVerIndex[2] = curTet->m_v[2];
		curFaceVerIndex[3] = curTet->m_v[3];
		
		if (!SearchFaceIn_AllTets(curFaceVerIndex[0],curFaceVerIndex[1],curFaceVerIndex[2],itet))
		{
			curTet->m_triType[3] = FACETYPE_BDY;
		}
		else
		{
			curTet->m_triType[3] = FACETYPE_INN;
		}
		
		if (!SearchFaceIn_AllTets(curFaceVerIndex[0],curFaceVerIndex[1],curFaceVerIndex[3],itet))
		{
			curTet->m_triType[2] = FACETYPE_BDY;
		}
		else
		{
			curTet->m_triType[2] = FACETYPE_INN;
		}
		
		if (!SearchFaceIn_AllTets(curFaceVerIndex[0],curFaceVerIndex[2],curFaceVerIndex[3],itet))
		{
			curTet->m_triType[1] = FACETYPE_BDY;
		}
		else
		{
			curTet->m_triType[1] = FACETYPE_INN;
		}
		
		if (!SearchFaceIn_AllTets(curFaceVerIndex[1],curFaceVerIndex[2],curFaceVerIndex[3],itet))
		{
			curTet->m_triType[0] = FACETYPE_BDY;
		}
		else
		{
			curTet->m_triType[0] = FACETYPE_INN;
		}
	}
	return 0;
}

bool SearchFaceIn_AllTets(const int verIndex1,const int verIndex2,const int verIndex3,const int skipTetIndex)
{
	if (verList.size() == 0 || tetMESH.size() == 0)
		return false;

	tetElement* curTet;
	int tetVerIndex[4];
	int identicalTimes;
	int k = 0;
	for (int itet = 0; itet < m_numTet; itet ++)
	{
		identicalTimes = 0;
		curTet = &(tetMESH[itet]);
		tetVerIndex[0] = curTet->m_v[0];
		tetVerIndex[1] = curTet->m_v[1];
		tetVerIndex[2] = curTet->m_v[2];
		tetVerIndex[3] = curTet->m_v[3];

		for (k = 0; k < 4; k ++)
		{
			if (tetVerIndex[k] == verIndex1 || tetVerIndex[k] == verIndex2 || tetVerIndex[k] == verIndex3)
				identicalTimes ++;
		}

		if (identicalTimes == 3 && itet != skipTetIndex)
			return true;
	}

	return false;
}

int GenerateTetFaces(void)
{
	std::vector<triElement> facesVector;	
	int ivTri[3],ivOther;
	int possibleTriIndex;// if one or more face of current tet has already, reocrd the abs index of the face in the tetFACES array
	
	int itet = 0;
	for (itet = 0; itet < m_numTet; itet ++)
	{
		((itet % (m_numTet / 10)) == 0) ? (printf("processed %d tet, altogether %d.\n",itet,m_numTet)) : (0);
		ivTri[0] = 0; ivTri[1] = 1; ivTri[2] = 2; ivOther = 3;
		possibleTriIndex = SearchtetFACES(itet,ivTri,ivOther,facesVector);// to process the first face of current tet
		
		ivTri[0] = 0; ivTri[1] = 1; ivTri[2] = 3; ivOther = 2;
		possibleTriIndex = SearchtetFACES(itet,ivTri,ivOther,facesVector);// to process the first face of current tet
		
		ivTri[0] = 0; ivTri[1] = 2; ivTri[2] = 3; ivOther = 1;
		possibleTriIndex = SearchtetFACES(itet,ivTri,ivOther,facesVector);// to process the first face of current tet
		
		ivTri[0] = 1; ivTri[1] = 2; ivTri[2] = 3; ivOther = 0;
		possibleTriIndex = SearchtetFACES(itet,ivTri,ivOther,facesVector);// to process the first face of current tet
	}
	
	m_numTetFACES = facesVector.size();
	
	if (tetFACES.size() > 0)
		tetFACES.clear();
	
	triElement newFACE;
	try {
		tetFACES.assign(m_numTetFACES,newFACE);
	}
	catch(...)
	{
		perror("Error when allocate memory for tetFACES");
		return -1;
	}
	
	int iface = 0;
	for (iface = 0; iface < m_numTetFACES; iface ++)
	{
		tetFACES[iface] = facesVector[iface];
	}	
	
	facesVector.clear();
	return 0;
}

int SearchtetFACES(int curTetIndex,int ivTri[3],int ivOther,std::vector<triElement> &faceVector)
{
	if (tetMESH.size() == 0)
		return -1;
	
	tetElement* curTet = &(tetMESH[curTetIndex]);
	int v1 = curTet->m_v[ivTri[0]];
	int v2 = curTet->m_v[ivTri[1]];
	int v3 = curTet->m_v[ivTri[2]];
	
	triElement possibleTri;
	triElement* foundTri;// record the found face
	int numCurFaces = faceVector.size();
	
	if (numCurFaces == 0)// initial the first element in faceVector
	{
		// update the face information
		possibleTri.m_v[0] = v1;
		possibleTri.m_v[1] = v2;
		possibleTri.m_v[2] = v3;
		
		possibleTri.m_tet[0] = curTetIndex;
		possibleTri.m_type = curTet->m_triType[ivOther];	
		possibleTri.m_ssb[0] = ivOther;	
		ComputeTriNormal3D(possibleTri);
		faceVector.push_back(possibleTri);
		
		// update the tet information
		curTet->m_faces[ivOther] = numCurFaces;
		
		return 0;
	}
	
	
	triElement* curFace;
	int vIdx1,vIdx2,vIdx3;	
	int iface = 0;
	for (iface = 0; iface < numCurFaces; iface ++)
	{
		curFace = &(faceVector[iface]);
		vIdx1 = curFace->m_v[0]; vIdx2 = curFace->m_v[1]; vIdx3 = curFace->m_v[2];
		
		if ((v1 == vIdx1 && v2 == vIdx2 && v3 == vIdx3) ||
			(v1 == vIdx1 && v2 == vIdx3 && v3 == vIdx2) ||
			(v1 == vIdx2 && v2 == vIdx1 && v3 == vIdx3) ||
			(v1 == vIdx2 && v2 == vIdx3 && v3 == vIdx1) ||
			(v1 == vIdx3 && v2 == vIdx1 && v3 == vIdx2) ||
			(v1 == vIdx3 && v2 == vIdx2 && v3 == vIdx1))
			
			break;// since there will not be a face which is adjacent to more than two tets, here we can break
		// iface stores the index of the searched face.
	}	
	
	if (iface == numCurFaces)// the v1v2v3 is not in the current faceVector, so we create a new triangle
	{
		// update the face information
		possibleTri.m_v[0] = v1;
		possibleTri.m_v[1] = v2;
		possibleTri.m_v[2] = v3;
		
		possibleTri.m_tet[0] = curTetIndex;
		possibleTri.m_type = curTet->m_triType[ivOther];
		possibleTri.m_ssb[0] = ivOther;
		ComputeTriNormal3D(possibleTri);		
		faceVector.push_back(possibleTri);
		
		// update the tet information
		curTet->m_faces[ivOther] = numCurFaces;
	}
	else// means the v1v2v3 face is already in the faceVector
	{
		// update the face info. some information should exist beforehand
		foundTri = &(faceVector[iface]);// don't know if this pointer will change the original data just like other pointer
		foundTri->m_tet[1] = curTetIndex;
		foundTri->m_ssb[1] = ivOther;
		curTet->m_faces[ivOther] = iface;
		curTet->m_tets[ivOther] = foundTri->m_tet[0];
		tetMESH[foundTri->m_tet[0]].m_tets[foundTri->m_ssb[0]] = curTetIndex;
	}
	
	return iface;
}

int ComputeTriNormal3D(triElement& tri)
{
	if (tri.m_type == -1)
		return -1;
	DBL3DVECT V1,V2,V3;
	DBL3DVECT DV1,DV2;
	V1 = VerElement2DBLVECT(verList[tri.m_v[0]]);
	V2 = VerElement2DBLVECT(verList[tri.m_v[1]]);
	V3 = VerElement2DBLVECT(verList[tri.m_v[2]]);
	
	DV1 = V2 - V1;
	DV2 = V3 - V1;
	
	DBL3DVECT crossProduct = DV1 / DV2;
	
	double vecLength = crossProduct.Len();
	
	if (fabs(vecLength) < 1e-10)
	{
		perror("There is degenerate triangle!");
		return -1;
	}
	else
	{
		tri.m_nx = crossProduct.x / vecLength;
		tri.m_ny = crossProduct.y / vecLength;
		tri.m_nz = crossProduct.z / vecLength;
	}
	
	// to make the normal point to the inner of the object
	DBL3DVECT triCenter = (V1 + V2 + V3) / 3;
	
	DBL3DVECT oppositePoint;
	oppositePoint = VerElement2DBLVECT(verList[tetMESH[tri.m_tet[0]].m_v[tri.m_ssb[0]]]);
	
	DBL3DVECT dirVec;
	dirVec = oppositePoint - triCenter;
	
	vecLength = dirVec.Len();
	if (vecLength < 1e-10)
	{
		perror("There is degenerate tet element");
		return -1;
	}
	
	dirVec = dirVec / vecLength;
	
	if (dirVec.x * tri.m_nx + dirVec.y * tri.m_ny + dirVec.z * tri.m_nz < 0)
	{
		tri.m_nx = -tri.m_nx;
		tri.m_ny = -tri.m_ny;
		tri.m_nz = -tri.m_nz;
	}
	
	return 0;
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

int InitBasicData(void)
{	
	m_numVer = m_numTet = m_numTetFACES /*= m_numTri*/ = 0;
	
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
	if (strcmp(fileFormat,"TET") != 0)
	{
		perror("File format is not TET.");
		return -1;
	}	
	
	fscanf(fp,"%d%d",&m_numVer,&m_numTet);	
	
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
	newVer->m_x = x;
	newVer->m_y = y;
	newVer->m_z = z;
	
	verList[0] = newVer;
	int iver = 0;
	for (iver = 1; iver < m_numVer; iver ++)
	{
 		newVer = new verElement[1];
		fscanf(fp,"%lf%lf%lf",&x,&y,&z);
		newVer->m_x = x; newVer->m_y = y; newVer->m_z = z;
		verList[iver] = newVer;
	}

	int itet = 0;
	for (itet = 0; itet < m_numTet; itet ++)
	{
		fscanf(fp,"%d %d%d%d%d %d%d%d%d",&meshType,\
			&(tetMESH[itet].m_v[0]),&(tetMESH[itet].m_v[1]),&(tetMESH[itet].m_v[2]),&(tetMESH[itet].m_v[3]),\
			&(tetMESH[itet].m_triType[0]),&(tetMESH[itet].m_triType[1]),&(tetMESH[itet].m_triType[2]),&(tetMESH[itet].m_triType[3]));		
	}

	fclose(fp);

	ComputeTetInfo_forVer();
		
	ComputeBdyInfo_forTet();

	CorrectTetInfo();
	ComputeAllFaceNormal();
	
	ComputeBdyVerNormal();

	// statistic the number of BOUNDARY faces
	int faceNum = 0,kFace = 0;
	const int numTetFaces = tetFACES.size();
	for (kFace = 0; kFace < numTetFaces; kFace ++)
	{
		if (tetFACES[kFace].m_tet[1] == -1)
			faceNum ++;
	}
		
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
	int kVer = 0;
	for (kVer = 0; kVer < m_numVer; kVer ++)
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
	int itet = 0;
	for (itet = 0; itet < m_numTet; itet ++)
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
	int iface = 0;
	for (iface = 0; iface < m_numTetFACES; iface ++)
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
	int iFace = 0;
	for (iFace = 0; iFace < m_numTetFACES; iFace ++)
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

double Area(const DBL3DVECT& P1,const DBL3DVECT& P2,const DBL3DVECT& P3)
{// the area of the triangle formed by P1,P2,P3
	return ((P2 - P1) / (P3 - P1)).Len() / 2;
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

	if (tetFACES.size() > 0)
		tetFACES.clear();

	return 0;
}
