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


#include "stdio.h"
#include "DataStruct.h"

DBL3DVECT::DBL3DVECT()
{
	x = y = z = 0;
	m_normal[0] = m_normal[1] = m_normal[2] = 0;
}

DBL3DVECT::DBL3DVECT(const double rx,const double ry,const double rz)
{
	x = rx; y = ry; z = rz;
}
DBL3DVECT::~DBL3DVECT()
{

}

bool DBL3DVECT::operator==(const DBL3DVECT& right)const
{
	if (fabs(x - right.x) + fabs(y - right.y) + fabs(z - right.z) < 1e-10)
		return true;
	else
		return false;
}

DBL3DVECT DBL3DVECT::operator+ (const DBL3DVECT& right)const
{
	DBL3DVECT sum;
	sum.x = x + right.x;
	sum.y = y + right.y;
	sum.z = z + right.z;
	return sum;
}

void DBL3DVECT::operator +=(const DBL3DVECT& right)
{
	this->x += right.x;
	this->y += right.y;
	this->z += right.z;
}

DBL3DVECT DBL3DVECT::operator - (const DBL3DVECT& right)const
{
	DBL3DVECT sub;
	sub.x = x - right.x;
	sub.y = y - right.y;
	sub.z = z - right.z;
	return sub;
}

DBL3DVECT DBL3DVECT::operator *(double factor)const
{
	DBL3DVECT pro;
	pro.x = x * factor;
	pro.y = y * factor;
	pro.z = z * factor;
	return pro;
}

double DBL3DVECT::operator *(const DBL3DVECT& right)const
{
	return this->x * right.x + this->y * right.y + this->z * right.z;
}

DBL3DVECT DBL3DVECT::operator /(double factor)const
{
	if (fabs(factor) < 1e-10)
	{
		perror("DBL3DVECT divided by 0.");
		return *this;
	}

	DBL3DVECT div;
	div.x = x / factor;
	div.y = y / factor;
	div.z = z / factor;
	return div;
}

void DBL3DVECT::operator /=(double factor)
{
	this->x /= factor;
	this->y /= factor;
	this->z /= factor;
}

DBL3DVECT DBL3DVECT::operator /(const DBL3DVECT& right)const
{
	DBL3DVECT crossProduct;
	crossProduct.x = this->y * right.z - this->z * right.y;
	crossProduct.y = -this->x * right.z + this->z * right.x;
	crossProduct.z = this->x * right.y - this->y * right.x;

	return crossProduct;
}

DBL3DVECT DBL3DVECT::operator -(void)const
{
	DBL3DVECT tmp;
	tmp.x = -this->x;
	tmp.y = -this->y;
	tmp.z = -this->z;
	return tmp;
}

double DBL3DVECT::Len()const
{
	return sqrt(x * x + y * y + z * z);
}

//////////////////////////////////////////////////////////////////////////
TetMesh::~TetMesh()
{

}

//////////////////////////////////////////////////////////////////////////
void edgeElement::Init()
{
	m_v[0] = m_v[1] = -1;
	m_type = EDGETYPE_INN;
	m_numTets = 0;
	m_Tets = NULL;
	m_ssb[0] = m_ssb[1] = NULL;
	m_ssbNB[0] = m_ssbNB[1] = NULL;
}

edgeElement::edgeElement()
{
	Init();
}

void edgeElement::Release()
{
	RELEASEMEM(m_Tets);
	RELEASEMEM(m_ssb[0]);
	RELEASEMEM(m_ssb[1]);
	RELEASEMEM(m_ssbNB[0]);
	RELEASEMEM(m_ssbNB[1]);
}

edgeElement::~edgeElement()
{
	Release();
}


//////////////////////////////////////////////////////////////////////////
verElement::verElement()
{
	m_x = m_y = m_z = 0;
	m_type = -1;
	m_numTet = 0;
	m_tets = NULL;
	m_subscriptIn_tet = NULL;
	m_subscriptIn_tet_2nd = NULL;
	m_subscriptIn_tet_3rd = NULL;
	m_subscriptIn_tet_4th = NULL;

	m_normal[0] = m_normal[1] = m_normal[2] = 0.0f;
	m_displace = 0;
}

verElement::verElement(const double x,const double y,const double z)
{
	verElement();
	m_x = x;
	m_y = y;
	m_z = z;
}

verElement::~verElement()
{
	if (m_tets != NULL)
		delete[] m_tets;
	m_tets = NULL;

	if (m_subscriptIn_tet != NULL)
		delete[] m_subscriptIn_tet;
	m_subscriptIn_tet = NULL; 

	if (m_subscriptIn_tet_2nd != NULL)
		delete[] m_subscriptIn_tet_2nd;
	m_subscriptIn_tet_2nd = NULL;

	if (m_subscriptIn_tet_3rd != NULL)
		delete[] m_subscriptIn_tet_3rd;
	m_subscriptIn_tet_3rd = NULL;

	if (m_subscriptIn_tet_4th != NULL)
		delete[] m_subscriptIn_tet_4th;
	m_subscriptIn_tet_4th = NULL;
}

verElement& verElement::operator =(const verElement& rightHand)
{
	this->m_x = rightHand.m_x;
	this->m_y = rightHand.m_y;
	this->m_z = rightHand.m_z;

	return (*this);
}

DBL3DVECT verElement::operator +(const verElement& rightHand)const 
{
	DBL3DVECT SUM;
	SUM.x = this->m_x + rightHand.m_x;
	SUM.y = this->m_y + rightHand.m_y;
	SUM.z = this->m_z + rightHand.m_z;

	return SUM;
}

DBL3DVECT verElement::operator -(const verElement& rightHand)const
{
	DBL3DVECT SUM;
	SUM.x = this->m_x - rightHand.m_x;
	SUM.y = this->m_y - rightHand.m_y;
	SUM.z = this->m_z - rightHand.m_z;

	return SUM;
}

DBL3DVECT verElement::operator *(double factor)const
{
	DBL3DVECT PRODUCT;
	PRODUCT.x = this->m_x * factor;
	PRODUCT.y = this->m_y * factor;
	PRODUCT.z = this->m_z * factor;

	return PRODUCT;
}

DBL3DVECT verElement::operator /(double factor)const
{
	DBL3DVECT PRODUCT;
	PRODUCT.x = this->m_x / factor;
	PRODUCT.y = this->m_y / factor;
	PRODUCT.z = this->m_z / factor;

	return PRODUCT;
}

//////////////////////////////////////////////////////////////////////////
triElement::triElement()
{
	m_v[0] = m_v[1] = m_v[2] = -1;
	m_type = -1;
	m_tet[0] = m_tet[1] = -1;
	m_ssb[0] = m_ssb[1] = -1;
	m_nx = m_ny = m_nz = 0;
}

triElement::~triElement()
{
	
}

//////////////////////////////////////////////////////////////////////////
tetElement::tetElement()
{
	m_v[0] = m_v[1] = m_v[2] = m_v[3] = -1;
	m_type = -1;
	m_faces[0] = m_faces[1] = m_faces[2] = m_faces[3] = -1;
	m_tets[0] = m_tets[1] = m_tets[2] = m_tets[3] = -1;
	m_Edges[0] = m_Edges[1] = m_Edges[2] = m_Edges[3] = m_Edges[4] = m_Edges[5] = -1;
	m_triType[0] = m_triType[1] = m_triType[2] = m_triType[3] = -1;
	
	m_minDiAngle = m_maxDiAngle = -1;
	m_bAngleComputed = false;
}

tetElement::~tetElement()
{
}