%{
#include "WrapSurfaceMesh.h"
%}

%feature("python:slot", "tp_str", functype="reprfunc") WrappedSimplicialComplex::as_string;
%include "WrapSurfaceMesh.h"
