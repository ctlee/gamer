# Geometry-preserving Adaptive MeshER
[![DOI](https://zenodo.org/badge/122682242.svg)](https://zenodo.org/badge/latestdoi/122682242)

GAMer is a surface mesh improvement library developed to condition surface meshes derived from noisy biological imaging data.
Using Tetgen, GAMer can generate tetrahedral meshes suitable for finite elements simulations of reaction-diffusion systems among others.
Currently this version of GAMer is available only as a library, Blender addon, and python module.
In the future, we will reintroduce CLI programs for mesh processing.

## Build Status
[![Master Build Status](https://travis-ci.org/ctlee/gamer.svg?branch=master)](https://travis-ci.org/ctlee/gamer)
[![Developement Build Status](https://travis-ci.org/ctlee/gamer.svg?branch=development)](https://travis-ci.org/ctlee/gamer)

## Installing
Prebuilt binaries of the GAMer Blender-addon are available under [releases](https://github.com/ctlee/gamer/releases).
Download the corresponding `.zip` for your platform and follow Blender's instructions to [install from file](https://docs.blender.org/manual/fi/dev/preferences/addons.html#header).

### Prerequisites
To build the GAMer library you will need access to a working C++ compiler with full C++14 support.
If you wish to use build and use the GAMer Python extensions, you will also need [SWIG > 3.0](http://www.swig.org/), access to a python intepreter, and the corresponding python shared library (`python.so` or `python.dylib`).
In order to use the GAMer Blender addon you should also have a working installation of Blender.

### Building
The following instructions are to build the base GAMer library.
If you wish to additionally compile the Blender GAMer addon, GAMer documentation, or other features please refer to the Additional Options section prior to building.

First, download a copy of the source from [releases](https://github.com/ctlee/gamer/releases) or clone the master branch.\
```bash
git clone https://github.com/ctlee/gamer.git
cd gamer
git submodule init
git submodule update
```

Linux and Mac:
```bash
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr -DBUILD_TESTS=on -DCMAKE_BUILD_TYPE=Release ..
make
```

For Windows, we support building using Microsoft Visual Studio (MSVS) through the use of CMake generators:
```bash
mkdir build64
cd build64
cmake -DBUILD_BLENDER=TRUE -G "Visual Studio 15 2017 Win64" ..
cmake --build . --config Release
```
If you get an `ImportError: DLL load failed` you are likely linking a different python library version than Blender's bundled python.
We recommend using Anaconda to obtain a python version matching Blender.

### Additional Options
To enable these options append the flags to your initial CMake function call
(e.g., `cmake -DBUILD_BLENDER=on ..`).
These can be used in addition to the standard [CMake flags](https://cmake.org/cmake/help/latest/manual/cmake.1.html).

* Build Python Extension: `-DBUILD_PYTHONEXT=ON`
* Build the Blender addon: `-DBUILD_BLENDER=ON`\
    This automatically builds the Python extension.
* Single Precision floating point numbers: `-DSINGLE=ON`
* Compile Tetgen Binary: `-DBUILD_TETGEN_BIN=ON`
* Test Cases: `-DBUILD_TESTS=ON`
    The test suite can be run using `make tests` or executing the output `./bin/objecttests`. Tests are organized using [GoogleTest](https://github.com/google/googletest) project.
* The GAMer documentation can be built if you have access to
[doxygen](http://www.doxygen.nl/) by running `make docs` with any build.

## Acknowledging your use of GAMer
Thanks for using GAMer! The developers would love to hear how you are using the tool. Please send us an email or post on GitHub letting us know. 

Please cite the above Zenodo DOI to acknowledge the software version and cite the following paper:\
[Lee, C. T.; Laughlin, J. G.; Angliviel de La Beaumelle, N.; Amaro, R.; McCammon, J. A.; Ramamoorthi, R.; Holst, M. J.; Rangamani, P. GAMer 2: A System for 3D Mesh Processing of Cellular Electron Micrographs. bioRxiv 2019, 534479.](https://www.biorxiv.org/content/10.1101/534479v1)

## Contributing
If you find a bug, unexpected behavior, or otherwise have some feature you would like to see please file an [Issue on GitHub](https://github.com/ctlee/gamer/issues).
To contribute to the development of GAMer please fork this repository and submit a pull request describing your additions.

## Authors
**[Christopher Lee](https://github.com/ctlee)**\
Department of Chemistry & Biochemistry\
University of California, San Diego

**John Moody**\
Department of Mathematics\
University of California, San Diego

### Contributors to GAMer
* Zeyun Yu (UCSD) and Yuhui Cheng(UCSD)\
Development of GAMer v1. To acknowledge your use of GAMer 1, please cite:\
[Yu, Z.; Holst, M. J.; Cheng, Y.; McCammon, J. A. Feature-Preserving Adaptive Mesh Generation for Molecular Shape Modeling and Simulation. J. Mol. Graph. Model. 2008, 26 (8), 1370–1380.](https://doi.org/10.1016/j.jmgm.2008.01.007)

* Tom Bartol (Salk Institute) and Johan Hake\
Development of Blender GAMer addon.

See also the list of [contributors](https://github.com/ctlee/gamer/contributors) who participated in this project.

## External libraries bundled with GAMer
* GAMer uses [Tetgen](http://wias-berlin.de/software/tetgen/) to generate
tetrahedralizations.

* GAMer uses [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) to
compute eigenvalues/eigenvectors of the Local Structure Tensor.

* GAMer uses [casc](https://github.com/ctlee/casc) as the underlying simplicial
complex data structure.

* GAMer uses [GoogleTest](https://github.com/google/googletest) to handle testing.

* Mesh checks in the GAMer Blender addon are inspired or borrowed from 3D Print Toolbox by Campbell Barton and Meshalyzer from CellBlender.

* [Triangle](https://www.cs.cmu.edu/~quake/triangle.html) is also bundled with GAMer but not currently used.
