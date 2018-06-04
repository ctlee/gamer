
%module (package="gamer") gamer

%{
#include "gamer.h"
%}

%include <std_string.i>

%include "typemaps.i"
%include "exceptions.i"
%include "std_unique_ptr.i"

%include "Vertex.i"
%include "SurfaceMesh.i"