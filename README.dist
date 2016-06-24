Release snapshot of GAMer
=========================

GAMer is a surface mesh improvement library included in the FEtk 
software umbrella. It depends on Maloc a Minimal Abstraction Layer 
for Object-oriented C, which all modules in FEtk depends on. Maloc 
is provided in this snapshot.


Installation
------------
This installation instructions should work on any UNIX based platform.

  # Only needed during install
  export PREFIX=installation/path
  export FETK_INCLUDE=$PREFIX/include
  export FETK_LIBRARY=$PREFIX/lib

  cd maloc
  configure --prefix=$PREFIX
  make
  make install

  cd gamer
  configure --prefix=$PREFIX
  make
  make install

Installation of PyGAMer
~~~~~~~~~~~~~~~~~~~~~~~

  cd gamer/swig
  configure --prefix=$PREFIX
  make
  make install

To compile the Python extension module of GAMer you need to have SWIG, 
the Python header files and NumPy installed.
 
Make sure to set LD_LIBRARY_PATH or (DYLD_LIBRARY_PATH on Mac) and 
PYTHONPATH before you run any programs or try using PyGAMer

  export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
  export PYTHONPATH=$PREFIX/lib/python2.6/site-packages:$PYTHONPATH

On Mac:

  export DYLD_LIBRARY_PATH=$PREFIX/lib:$DYLD_LIBRARY_PATH
  export PYTHONPATH=$PREFIX/lib/python2.6/site-packages:$PYTHONPATH

Installation of GAMer applications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  cd gamer/tools/{ImproveSurfMesh, MolecularMesh, GenerateMesh}
  configure --prefix=$PREFIX
  make
  make install

PyGAMer
-------
PyGAMer is a Python wrapper of the core GAMer library. To test the
main functionality you can run the test script in gamer/swig/test.

  cd gamer/swig/test
  python gamer_test.py

Blender
-------
A blender plug-in for GAMer is provided in:

  gamer/tools/blender

In addition to following the instructions in the README file in the 
blender directory you also need to have a functional PyGAMer installation.

Contact
-------
Please report any problems to:

  Johan Hake <hake.dev@gmail.com>

San Diego 2/16/11
