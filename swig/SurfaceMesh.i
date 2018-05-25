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

%include "std_unique_ptr.i"


wrap_unique_ptr(SMUniquePtr, SurfaceMesh);

%include "SurfaceMesh.h"

%inline %{
class StopIterator {};

template <typename IT, typename T>
class IteratorWrapper {
public:
    IteratorWrapper(IT _begin, IT _end): curr(_begin), end(_end) {}
	
	IteratorWrapper* __iter__()
    {
      return this;
    }

    // Iterator function for python 2.x
    T& next(){
    	if(curr != end)
    		return *curr++;
    	else
    		throw StopIterator();
    }

    // Iterator function for python 3.x
    T& __next__(){
    	if(curr != end)
    		return *curr++;
    	else
    		throw StopIterator();
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

class Face {};

// %extend IteratorWrapper<SMVDataIterator, Vertex>
// {
//   Vertex& next()
//   {
//     if ($self->curr != $self->end)
//     {
//       // dereference the iterator and return reference to the object,
//       // after that it increments the iterator
//       return *$self->curr++;
//     }
//     throw StopIterator();
//   }
// }


%template(VIT) IteratorWrapper<SMVDataIterator, Vertex>;
%template(FIT) IteratorWrapper<SMFDataIterator, Face>;

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

	%pythoncode %{
		def vertices(self):
			for v in self.getVertexIT():
				yield v
	%}
	}
};

%rename(printMesh) print(const SurfaceMesh &mesh);
void print(const SurfaceMesh& mesh);

// import _gamer as g
// mesh = g.SurfaceMesh()
// mesh.insertVertex(g.VertexKey([1]), g.Vertex(1,2,3))