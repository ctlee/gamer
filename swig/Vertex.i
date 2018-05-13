%{
#include "Vertex.h"
%}
%ignore Vertex::operator[](std::size_t index);
%ignore Vertex::operator[](std::size_t index) const;
%ignore Vertex::operator=(const Vertex& v); //created a separate copy function if necessary
%ignore Vertex(const Vertex& x);
%ignore Vertex(const Vertex&& x);


%rename(vector) operator Vector;
%rename(assign) operator=(const Vertex& v);
%rename(boolequal) operator==(const Vertex& rhs);
%rename(add_vector) operator+(const Vertex& A, const Vector& B);
%rename(add_vertices) operator+(const Vertex& A, const Vertex& B);
%rename(subtract_vector) operator-(const Vertex& A, const Vector& B);
%rename(subtract_vertices) operator-(const Vertex& A, const Vertex& B);
%rename(scalar_multiply) operator*(double x, const Vertex& A);
%rename(scalar_divide) operator/(const Vertex& A, double x);
%rename(getitem) operator[](std::size_t index);
%rename(copy) operator=(const Vertex& v);

// %feature("python:slot", "tp_str", functype="reprfunc") Vertex::to_string;
%ignore operator<<(std::ostream& output, const Vertex& v);
%include "Vertex.h"
%feature("python:slot", "tp_repr", functype="reprfunc") Vertex::to_string;

// %extend Vertex{
// 	std::string to_string(){
// 		return $self->to_string();
// 	}	
// }


class Vector{
public:	
	Vector();
	~Vector();

	std::string as_string();
};

%extend Vertex{
	const double& getitem(std::size_t index) const {
		return (*($self))[index];
	}
	//void copy(Vertex& v) {
	//	Vertex vnew;
	//	vnew.position = v.position;
	//	vnew.marker = v.marker;
	//	vnew.selected = v.selected;
	//}
}