%{
#include "WrapSurfaceMesh.h"
#include "SurfaceMesh.h"
%}

%include <std_array.i>
%include <carrays.i>
namespace std {
	%template(VertexKey) array<int,1>;
	%template(EdgeKey) array<int,2>;
	%template(FaceKey) array<int, 3>;
	//%template(Edge) array<Vertex, 2>;
}

%include "Vertex.i"

%rename(_print_) print(const SurfaceMesh &mesh);

%feature("python:slot", "tp_str", functype="reprfunc") WrappedSimplicialComplex::as_string;
%include "WrapSurfaceMesh.h"

class SurfaceMesh{

public:	
	SurfaceMesh();
	~SurfaceMesh();


	%extend {
		void insertVertex(std::array<int, 1> &s, const Vertex &data) {
			$self->insert(s, data);
		}
	}

	%extend {
		void insertEdge(std::array<int, 2> &s) {
			$self->insert(s);
		}
	}


	%extend {
		void insertFace(std::array<int, 3> &s, const Face &data) {
			$self->insert(s, data);
		}
	}

};

void print(const SurfaceMesh& mesh);