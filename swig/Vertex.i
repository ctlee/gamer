%{
#include "Vertex.h"
%}

%rename(add_vertices) operator+(const Vertex& A, const Vertex& B);

%include "Vertex.h"


class Vector{
public:	
	Vector();
	~Vector();

	std::string as_string();
};