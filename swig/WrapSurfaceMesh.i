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
}

%include "Vertex.i"

%rename(_print_) print(const SurfaceMesh &mesh);

%feature("python:slot", "tp_str", functype="reprfunc") WrappedSimplicialComplex::as_string;
%include "WrapSurfaceMesh.h"

class SurfaceMesh{


//%rename (insert_vertex) insert<1>(const std::array<int, 1> &s, const Vertex &data);


public:	
	SurfaceMesh();
	~SurfaceMesh();

	%extend {
		void insertVertex(std::array<int, 1> &s, const Vertex &data) {
			$self->insert(s, data);
		}
	}

	//void insert<1>(const std::array<int, 1> &s, const Vertex &data);
};

void print(const SurfaceMesh& mesh);