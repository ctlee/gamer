Release snapshot of GAMer
=========================

GAMer is a surface mesh improvement library included in the FEtk 
software umbrella. It depends on Maloc a Minimal Abstraction Layer 
for Object-oriented C, which all modules in FEtk depends on. Maloc 
is provided in this snapshot.

Quick Installation
------------------

```bash
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr -DCMAKE_BUILD_TYPE=Release ..
make all
make install
```

Blender
-------
A blender plug-in for GAMer is provided in:

  gamer/tools/blender
  
To build it append to your intial CMake configuration.

```bash
-DBUILD_BLENDER=ON
```

In addition to following the instructions in the README file in the 
blender directory you also need to have a functional PyGAMer installation.

Full Installation
-----------------

Tests
-----
```
-DBUILD_TESTS=ON
```

Contact
-------
Please report any problems to:

    Christopher Lee <ctlee@ucsd.edu>
    John Moody <jbmoody@ucsd.edu>

