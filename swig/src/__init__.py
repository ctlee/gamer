""" Python interface of GAMer
"""

__author__ = "Johan Hake (hake.dev@gmail.com)"
__copyright__ = "Copyright (C) 2010 Johan Hake"
__date__ = "2010-08-06 -- 2012-04-18"
__license__  = "GNU LGPL Version 3.1"

# Import the wrapped SWIG wrapped version of gamer and the upy gui
import gamer.cgamer

# Expose 2 classes from the C interface of gamer
SurfaceMesh = gamer.cgamer.SurfaceMesh
GemMesh = gamer.cgamer.GemMesh

# PDB read functions returning a SurfaceMesh
read_PDB_molsurf = gamer.cgamer.SurfaceMesh_readPDB_molsurf
read_PDB_gauss = gamer.cgamer.read_PDB_gauss
