%{
#include "SurfaceMesh.h"
#include "libraries/casc/casc"

%}

%include <std_array.i>
%include <carrays.i>
namespace std {
	%template(VertexKey) array<int,1>;
	%template(EdgeKey) array<int,2>;
	%template(FaceKey) array<int, 3>;
}

%include "Vertex.i"

%inline %{
class StopIterator {};

template <typename IT, typename T>
class IteratorWrapper {
public:
    Iterator(IT _cur, IT _end) : cur(_cur), end(_end) {}
    Iterator* __iter__()
    {
      return this;
    }
    IT cur;
    IT end;
  };
%}

%include "exception.i"
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

%extend Iterator
{
  T& next()
  {
    if ($self->cur != $self->end)
    {
      // dereference the iterator and return reference to the object,
      // after that it increments the iterator
      return *$self->cur++;
    }
    throw StopIterator();
  }
}

// %template(VertexIterator) Iterator<VertexIT, Vertex>;

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

		// VertexIterator getVertexIT(){
		// 	return std::static_cast<VertexIterator>($self->get_level<1>());
		// }

	// %pythoncode%{
	// 	def vertices(self):
	// 		iterable = self.getVertexIT();	
	// 		for i in range(self.sizeVertices()):
	// 			yield next(iterable)
	// %}
	}
};

// def vertices(self):
//     "Return an iterator over vertices"
//     for i in range(self.num_vertices):
//         yield _Vertex_getitem(self, i)

%rename(printMesh) print(const SurfaceMesh &mesh);
void print(const SurfaceMesh& mesh);

// import _gamer as g
// mesh = g.SurfaceMesh()
// mesh.insertVertex(g.VertexKey([1]), g.Vertex(1,2,3))