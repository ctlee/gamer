
%module (package="gamer") gamer

%{
#include "gamer.h"
%}

%include <std_string.i>
%include <std_vector.i>

%include "typemaps.i"
%include "exceptions.i"
%include "std_unique_ptr.i"

%include "Vertex.i"
%include "SurfaceMesh.i"

%{
struct index_map {
  std::vector<int> map;

  index_map(std::size_t size){
    map = std::vector<int>(size);
  }

  void updatemap(int marker, PyObject* incoming){
   if (PyList_Check(incoming)) {
      for(Py_ssize_t i = 0; i < PyList_Size(incoming); i++) {
        PyObject *value = PyList_GetItem(incoming, i);
        map[PyLong_AsSize_t(value)] = marker;
      }
    }
    else {
      throw std::logic_error("Passed PyObject pointer was not a list!");
    }
  }

  int operator[](std::size_t idx){
    return map[idx];
  }

  void printMap(){
    for(auto item : map){
      std::cout << item << ", ";
    }
    std::cout << std::endl;
  }
};
%}

struct index_map {
  index_map(std::size_t size);
  void updatemap(int marker, PyObject* incoming);
  void printMap();
  %extend{
    int __getitem__(std::size_t key){
      return $self->map[key];
    }
  }
};
