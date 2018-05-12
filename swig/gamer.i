
%module (package="gamer") gamer

%{
#include "gamer.h"
%}

%include <std_string.i>

%include "exceptions.i"
%include "Vertex.i"
%include "SurfaceMesh.i"