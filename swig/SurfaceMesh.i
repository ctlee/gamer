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
#include <iostream>
#include "SurfaceMesh.h"
#include "PDBReader.h"
%}

%inline %{
using VertexID_Iterator = SurfaceMesh::SimplexIDIterator<1>;
using VertexData_Iterator = SurfaceMesh::DataIterator<1>;

using FaceID_Iterator = SurfaceMesh::SimplexIDIterator<3>;
using FaceData_Iterator = SurfaceMesh::DataIterator<3>;

using VertexID = SurfaceMesh::SimplexID<1>;
using EdgeID = SurfaceMesh::SimplexID<2>;
using FaceID = SurfaceMesh::SimplexID<3>;
%}

%include <std_array.i>

%template(VertexKey) std::array<int,1>;
%template(EdgeKey) std::array<int,2>;
%template(FaceKey) std::array<int, 3>;

// Wrap the unique_ptr
wrap_unique_ptr(SMUniquePtr, SurfaceMesh);

struct Face
{
  Face() {}
  Face(int marker, bool selected) {} // Shadows full constructor
  int  marker;   /**< @brief Marker */
  bool selected; /**< @brief Selection flag */

  %pythoncode %{
    def __repr__(self):
      output = "Face(m:" + str(self.marker) \
          + ";sel:" + str(self.selected) + ")"
      return output
  %}
};

%include "SurfaceMesh.h"
%include "IteratorWrapper.i"

%template(VertexIT) IteratorWrapper<VertexData_Iterator, Vertex>;
%template(FaceIT) IteratorWrapper<FaceData_Iterator, Face>;

%template(FaceIDIT) IteratorWrapper<FaceID_Iterator, FaceID>;

class FaceID {
public:
  %extend {
    const char* __repr__() {
      std::ostringstream out;
      out << *$self;
      return out.str().c_str();
    }
  }
};


// Class to shadow complicated alias
class SurfaceMesh{
public:
  %extend {
    void addVertex(const Vertex &data){
      $self->add_vertex(data);
    }
    void insertVertex(std::array<int, 1> &s, const Vertex &data) {
      $self->insert<1>(s, data);
    }
    void insertEdge(std::array<int, 2> &s) {
      $self->insert<2>(s);
    }
    void insertFace(std::array<int, 3> &s, const Face &data) {
      $self->insert<3>(s, data);
    }
    void insertFace(std::array<int, 3> &s) {
      $self->insert<3>(s);
    }

    FaceID getFaceID(std::array<int, 3> &s) {
      return $self->get_simplex_up<3>(s);
    }

    Vertex& getVertex(std::array<int, 1> &s){
      return *($self->get_simplex_up<1>(s));
    }

    Face& getFace(std::array<int, 3> &s){
      return *($self->get_simplex_up<3>(s));
    }

    int sizeVertices(){
      return $self->size<1>();
    }
    int sizeEdges(){
      return $self->size<2>();
    }
    int sizeFaces(){
      return $self->size<3>();
    }

    std::array<int, 3> getFaceName(FaceID f){
      return $self->get_name(f);
    }

    IteratorWrapper<VertexData_Iterator, Vertex> getVertexIT(){
      auto it = $self->get_level<1>();
      return IteratorWrapper<VertexData_Iterator, Vertex>(it.begin(), it.end());
    }

    IteratorWrapper<FaceData_Iterator, Face> getFaceDataIT(){
      auto it = $self->get_level<3>();
      return IteratorWrapper<FaceData_Iterator, Face>(it.begin(), it.end());
    }

    IteratorWrapper<FaceID_Iterator, FaceID> getFaceIT(){
      auto it = $self->get_level_id<3>();
      return IteratorWrapper<FaceID_Iterator, FaceID>(it.begin(), it.end());
    }
  }

  %pythoncode %{
    def vertices(self):
      for v in self.getVertexIT():
        yield v

    def faces(self):
      for f in self.getFaceDataIT():
        yield f

    def faceIDs(self):
      for f in self.getFaceIT():
        yield f
  %}
};


%inline %{
SurfaceMesh* ReadOFF(const std::string &filename){
  return readOFF(filename).release();
}
SurfaceMesh* ReadOBJ(const std::string &filename){
  return readOBJ(filename).release();
}
SurfaceMesh* RefineMesh(const SurfaceMesh &mesh){
  return refineMesh(mesh).release();
}
SurfaceMesh* ReadPDB_molsurf(const std::string & filename){
  return readPDB_molsurf(filename).release();
}
SurfaceMesh* ReadPDB_gauss(const std::string& filename, float blobbyness, float isovalue){
  return readPDB_gauss(filename, blobbyness, isovalue).release();
}
%}

void writeOFF(const std::string &filename, const SurfaceMesh &mesh);
void writeOBJ(const std::string &filename, const SurfaceMesh &mesh);

bool smoothMesh(SurfaceMesh &mesh, int maxMinAngle, int minMaxAngle, int maxIter, bool preserveRidges);
void coarse(SurfaceMesh &mesh, double coarseRate, double flatRate, double denseWeight);
void scale(SurfaceMesh &mesh, double sx, double sy, double sz);
void translate(SurfaceMesh &mesh, double dx, double dy, double dz);
double getVolume(const SurfaceMesh &mesh);
double getArea(const SurfaceMesh &mesh);

%rename(printMesh) print(const SurfaceMesh &mesh);
void print(const SurfaceMesh &mesh);

void compute_orientation(SurfaceMesh& F); // Currently force return to void


%pythoncode %{
  def coarse_dense(self, mesh, rate=1.6, numiter=1):
    "Coarse flat areas"
    if not isinstance(rate, (int, float)) or rate <= 0:
      raise TypeError("expected a positive scalar for the 'rate' argument")
    if not isinstance(numiter, int) or numiter < 1:
      raise TypeError("expected a positive scalar for the 'numiter' argument")
    for i in range(numiter):
      coarse(self, rate, 0, 10) # coarse rate, flat rate, dense rate

  def coarse_flat(self, mesh, rate=0.016):
    "Coarse flat areas"
    coarse(self, mesh, rate, 1, 0) # coarse rate, flat rate, dense rate
%}