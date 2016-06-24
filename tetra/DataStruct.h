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

#ifndef DATASTURE_IN_TETSMOOTH
#define DATASTURE_IN_TETSMOOTH

#include "math.h"

class DBL3DVECT{
public:
	double x; double y; double z;
	double m_normal[3];
public:
	DBL3DVECT();
	DBL3DVECT(const double rx,const double ry,const double rz);
	~DBL3DVECT();
public:
	bool operator== (const DBL3DVECT& right)const;
	DBL3DVECT operator+ (const DBL3DVECT& right)const;
	void	  operator+=(const DBL3DVECT& right);
	DBL3DVECT operator- (const DBL3DVECT& right)const;
	DBL3DVECT operator* (double factor)const;
	double    operator* (const DBL3DVECT& right)const;
	DBL3DVECT operator/ (double factor)const;
	void	  operator/=(double factor);
	DBL3DVECT operator/(const DBL3DVECT& right)const;// the cross product of two vectors.
	DBL3DVECT operator-()const;

public:
	double Len(void)const;
};


class TetMesh{
public:
	class verElement;// vertex element
	class edgeElement;//edge element
	class tetElement;// tetrahedron element
	class triElement;// triangle element
			
public:
	TetMesh(){};
	~TetMesh();

	
	///////////////////////////////  Vertex Element ///////////////////////////////////////////
	class verElement{
	public:
		verElement();
		verElement(const double x,const double y,const double z);
		~verElement();
		
	public:
		double m_x,m_y,m_z;//coordinate of current vertex
		int m_type;//vertex type: boundary =0, 1 for layer1,  inner = 255 (-1 for default)
		//////////////////////////////////////////////////////////////////////////
		int m_numTet;// number tet
		int* m_tets;// index of all adjacent tets in tetMESH
		int* m_subscriptIn_tet;// the relative index of current vertex in the current tet. e.g. suppose m_numtet = 3,
		// then m_tets contains 3 tet and m_relativeI_tet contains 3 integers. and m_relativeI_tet[k]
		// is the subscript (0-3) of the current vertex in the tet[k]-th tetrahedra. i.e. verList[tetMESH[m_tet[k]].ver[m_relativeI_tet[k]]] is current vertex
		int* m_subscriptIn_tet_2nd;
		int* m_subscriptIn_tet_3rd;
		int* m_subscriptIn_tet_4th;
		double m_normal[3];
		double m_displace;// to denote the displace between current step and last step		
		
	public:
		verElement& operator=(const verElement& rightHand);
		DBL3DVECT  operator+(const verElement& rightHand)const;
		DBL3DVECT  operator-(const verElement& rightHand)const;
		DBL3DVECT  operator*(double factor)const;
		DBL3DVECT  operator/(double factor)const;
		inline double Mag(){return sqrt(m_x * m_x + m_y * m_y + m_z * m_z);};
	};
	
	/////////////////// Edge Element ///////////////////////////////////////////////////////
	class edgeElement{
	public:

		int m_v[2];// the abs index of the two end points of current edge
		edgeElement();
		~edgeElement();
		int m_type;// EDGETYPE_BDY or EDGETYPEINNER
		
	public:
		
		int m_numTets;// the number of tets which are adjacent to current edge
		int* m_Tets;// 
		int* m_ssb[2];// the relative index of the two end points in each tet
		int* m_ssbNB[2];// the relateive index of the two neighbors points in each tet
		
	public:
		void Init();
		void Release();
		
	};
	
	//////////////// Triangle Face Element //////////////////////////////////////////////////////////
	typedef struct{
		int V1;
		int V2;
		int V3;
	}SIMPLE_TRI;
	class triElement{
	public:
		triElement();
		~triElement();

		int m_v[3];// the index of the three vertex of current triangle
		int m_type;// boundary =0, inner 255,default -1
		int m_tet[2];// the two possible tets to which the current triangle may be adjacent (-1 means no tet)
		int m_ssb[2];// the subscripts of current triangle in the two possible tets 
		double m_nx,m_ny,m_nz;// the normal of current triangle. Science we show all the triangles in front and back face, 
		// the orientation of normal doesn't matter. 
		
	public:
		
	};
	
	////////////////// Tetrahedral Element ////////////////////////////////////////////////////////
	class tetElement{
	public:

		tetElement();
		~tetElement();
		
		int m_v[4]; // the index of the 4 vertex of current tet
		int m_triType[4];//the type of the 4 triangles of cur tet
		//	int m_bdytri[4];// the abs index in sufMESH of the 4 possible boundry triangles of cur tet. 
		
	public:// 
		int m_type;// boundary = 0, inner 255. And for future research, we can also code the tet by the its distance to the boundary in order to 
		// fast retrieves that kind of information.	
		
		
	public:
		int m_faces[4];// the abs index in tetFACES of the 4 faces.
		int m_tets[4];// the abs idx of the possible 4 adj tets
		int m_Edges[6];
		
	public:
		double m_minDiAngle;
		double m_maxDiAngle;// the min and max dihedral angle in this tet
		bool   m_bAngleComputed;
	};
};

typedef TetMesh::verElement verElement;
typedef TetMesh::edgeElement edgeElement;
typedef TetMesh::triElement triElement;
typedef TetMesh::tetElement tetElement;
typedef TetMesh::SIMPLE_TRI SIMPLE_TRI;

typedef struct {
	double x1;
	double y1;
	double z1;
	double x2;
	double y2;
	double z2;
	double x3;
	double y3;
	double z3;
}EIGENVECT;

#ifndef RELEASEMEM
#define RELEASEMEM(ptr) {\
	if (NULL != ptr){try{delete[] ptr;ptr = NULL;}catch(...){perror("RELEASE Error");}}\
}
#endif

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MATHPI
#define MATHPI (3.1415926535897)
#endif

#ifndef PTTYPE_DEF
#define PTTYPE_DEF
#define PTTYPE_ALL (-1)
#define PTTYPE_BDY 0
#define PTTYPE_INN 255
#define EDGETYPE_BDY PTTYPE_BDY
#define EDGETYPE_INN PTTYPE_INN
#define FACETYPE_BDY PTTYPE_BDY
#define FACETYPE_INN PTTYPE_INN
#define TETTYPE_BDY PTTYPE_BDY
#define TETTYPE_INN PTTYPE_INN
#endif

#endif