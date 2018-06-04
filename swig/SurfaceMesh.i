%{
#include <iostream>
#include "SurfaceMesh.h"

using SMVertIterator = SurfaceMesh::SimplexIDIterator<1>;
using SMVDataIterator = SurfaceMesh::DataIterator<1>;
using SMFDataIterator = SurfaceMesh::DataIterator<3>;
%}

%include <std_array.i>

%template(VertexKey) std::array<int,1>;
%template(EdgeKey) std::array<int,2>;
%template(FaceKey) std::array<int, 3>;

wrap_unique_ptr(SMUniquePtr, SurfaceMesh);
%rename(printMesh) print(const SurfaceMesh &mesh);

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
      std::cout << "Next()" << std::endl;
      T& obj = *curr;
      std::cout << "C++ object: " << obj <<  std::endl;
      std::cout << "C++ ptr: " << &obj << std::endl;
      ++curr;
      return obj;
    }
    else{
      std::cout << "Throwing StopIterator()" << std::endl;
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
%}
SurfaceMesh* ReadOFF(const std::string &filename);