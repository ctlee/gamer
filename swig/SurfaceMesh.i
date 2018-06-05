%{
#include <iostream>
#include "SurfaceMesh.h"

using SMVertIterator = SurfaceMesh::SimplexIDIterator<1>;
using SMVDataIterator = SurfaceMesh::DataIterator<1>;
using SMFDataIterator = SurfaceMesh::DataIterator<3>;
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
};

%include "SurfaceMesh.h"

%inline %{
class StopIterator {};

template <typename IT, typename T>
class IteratorWrapper {
public:
  IteratorWrapper(IT _begin, IT _end): curr(_begin), end(_end) {}

  IteratorWrapper* __iter__(){
    return this;
  }

  // Iterator function for Python 2.x
  T& next(){
    return __next__();
  }

  // Iterator function for Python 3.x
  T& __next__(){
    if(curr != end){
      T& obj = *curr;
      ++curr;
      return obj;
    }
    else{
      throw StopIterator();
    }
  }

  IT curr;
  IT end;
};
%}

%exception IteratorWrapper::next {
  try
  {
    $action // calls %extend function next() below
  }
  catch (StopIterator)
  {
    PyErr_SetString(PyExc_StopIteration, "End of iterator");
    return NULL;
  }
}

%exception IteratorWrapper::__next__ {
  try
  {
    $action // calls %extend function next() below
  }
  catch (StopIterator)
  {
    PyErr_SetString(PyExc_StopIteration, "End of iterator");
    return NULL;
  }
}

%template(VertexIT) IteratorWrapper<SMVDataIterator, Vertex>;
%template(FaceIT) IteratorWrapper<SMFDataIterator, Face>;

// Class to shadow complicated alias
class SurfaceMesh{
public:
  %extend {
    void insertVertex(std::array<int, 1> &s, const Vertex &data) {
      $self->insert(s, data);
    }
    void insertEdge(std::array<int, 2> &s) {
      $self->insert(s);
    }
    void insertFace(std::array<int, 3> &s, const Face &data) {
      $self->insert(s, data);
    }
    void insertFace(std::array<int, 3> &s) {
      $self->insert(s);
    }

    Vertex& getVertex(std::array<int, 1> &s){
      return *($self->get_simplex_up(s));
    }

    Face& getFace(std::array<int, 3> &s){
      return *($self->get_simplex_up(s));
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

    IteratorWrapper<SMVDataIterator, Vertex> getVertexIT(){
      auto it = $self->get_level<1>();
      return IteratorWrapper<SMVDataIterator, Vertex>(it.begin(), it.end());
    }

    IteratorWrapper<SMFDataIterator, Face> getFaceIT(){
      auto it = $self->get_level<3>();
      return IteratorWrapper<SMFDataIterator, Face>(it.begin(), it.end());
    }

  %pythoncode %{
    def vertices(self):
      for v in self.getVertexIT():
        yield v

    def faces(self):
      for f in self.getFaceIT():
        yield f
  %}
  }
};


%{
SurfaceMesh* ReadOFF(const std::string &filename){
  return readOFF(filename).release();
}
SurfaceMesh* ReadOBJ(const std::string &filename){
  return readOBJ(filename).release();
}
SurfaceMesh* RefineMesh(const SurfaceMesh &mesh){
  return refineMesh(mesh).release();
}
%}
SurfaceMesh* ReadOFF(const std::string &filename);
SurfaceMesh* ReadOBJ(const std::string &filename);
SurfaceMesh* RefineMesh(const SurfaceMesh &mesh);

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