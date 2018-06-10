// ***************************************************************************
// This file is part of the GAMer software.
// Copyright (C) 2016-2017
// by Christopher Lee, John Moody, Rommie Amaro, J. Andrew McCammon,
//    and Michael Holst

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

// ***************************************************************************

%{
#include "Vertex.h"
%}

%ignore Vertex::operator Vector; // Returns position
%ignore Vertex::operator[](std::size_t index);
%ignore Vertex::operator[](std::size_t index) const;
%ignore Vertex::operator=(const Vertex& v); //created a separate copy function if necessary
%ignore Vertex(const Vertex& x);
%ignore Vertex(const Vertex&& x);
%ignore Vertex::to_string() const; // Passing strings from c++ to python often segmentation faults...

// %ignore Vertex::operator=(const Vertex& v);
// %ignore Vertex::operator==(const Vertex& rhs);
%ignore operator+(const Vertex& A, const Vector& B);
%ignore operator+(const Vertex& A, const Vertex& B);
%ignore operator-(const Vertex& A, const Vector& B);
%ignore operator-(const Vertex& A, const Vertex& B);
%ignore operator*(double x, const Vertex& A);
%ignore operator/(const Vertex& A, double x);
%ignore operator=(const Vertex& v);

// %rename(assign) operator=(const Vertex& v);
// %rename(boolequal) operator==(const Vertex& rhs);
// %rename(add_vector) operator+(const Vertex& A, const Vector& B);
// %rename(add_vertices) operator+(const Vertex& A, const Vertex& B);
// %rename(subtract_vector) operator-(const Vertex& A, const Vector& B);
// %rename(subtract_vertices) operator-(const Vertex& A, const Vertex& B);
// %rename(scalar_multiply) operator*(double x, const Vertex& A);
// %rename(scalar_divide) operator/(const Vertex& A, double x);
// %rename(copy) operator=(const Vertex& v);

%ignore operator<<(std::ostream& output, const Vertex& v);
%include "Vertex.h"

%extend Vertex{
	const double& getitem(std::size_t index) const {
		return (*($self))[index];
	}

  const char* __repr__() {
    return $self->to_string().c_str();
  }
}