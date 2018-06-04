%{
#include "Vertex.h"
%}

%ignore Vertex::operator Vector; // Returns position
%ignore Vertex::operator[](std::size_t index);
%ignore Vertex::operator[](std::size_t index) const;
%ignore Vertex::operator=(const Vertex& v); //created a separate copy function if necessary
%ignore Vertex(const Vertex& x);
%ignore Vertex(const Vertex&& x);
%ignore Vertex::to_string() const; // Passing strings from c++ to python often segmentation faults...



%ignore operator=(const Vertex& v);
%ignore operator==(const Vertex& rhs);
%ignore operator+(const Vertex& A, const Vector& B);
%ignore operator+(const Vertex& A, const Vertex& B);
%ignore operator-(const Vertex& A, const Vector& B);
%ignore operator-(const Vertex& A, const Vertex& B);
%ignore operator*(double x, const Vertex& A);
%ignore operator/(const Vertex& A, double x);
%ignore operator=(const Vertex& v);

// %rename(assign) operator=(const Vertex& v);
// %rename(boolequal) operator==(const Vertex& rhs);
// %rename(add_vector) operator+(const Vertex& A, const Vector& B);
// %rename(add_vertices) operator+(const Vertex& A, const Vertex& B);
// %rename(subtract_vector) operator-(const Vertex& A, const Vector& B);
// %rename(subtract_vertices) operator-(const Vertex& A, const Vertex& B);
// %rename(scalar_multiply) operator*(double x, const Vertex& A);
// %rename(scalar_divide) operator/(const Vertex& A, double x);
// %rename(copy) operator=(const Vertex& v);

%ignore operator<<(std::ostream& output, const Vertex& v);
%include "Vertex.h"

// class Vector{
// public:
// 	Vector();
// 	~Vector();

// 	std::string as_string();
// };

%extend Vertex{
	const double& getitem(std::size_t index) const {
		return (*($self))[index];
	}

  %pythoncode %{
    def __repr__(self):
      output = "Vertex(x:" + str(self.getitem(0)) \
          + ";y:" + str(self.getitem(1)) \
          + ";z:" + str(self.getitem(2)) \
          + ";m:" + str(self.marker) + " sel:" + str(self.selected) + ")"
      return output
  %}
}