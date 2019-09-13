# Geometry-preserving Adaptive MeshER
[![DOI](https://zenodo.org/badge/122682242.svg)](https://zenodo.org/badge/latestdoi/122682242)
[![PyPI](https://img.shields.io/pypi/v/pygamer)](https://pypi.org/project/pygamer/)
[![Master Build Status](https://travis-ci.org/ctlee/gamer.svg?branch=master)](https://travis-ci.org/ctlee/gamer)
[![Build status](https://ci.appveyor.com/api/projects/status/urffu7062fnohidl/branch/master?svg=true)](https://ci.appveyor.com/project/ctlee/gamer)
[![Documentation Status](https://readthedocs.org/projects/gamer/badge/?version=latest)](https://gamer.readthedocs.io/en/latest/?badge=latest)

GAMer is a surface mesh improvement library developed to condition surface meshes derived from noisy biological imaging data.
Using Tetgen, GAMer can generate tetrahedral meshes suitable for finite elements simulations of reaction-diffusion systems among others.
GAMer has the following main features:

* Surface mesh improvement and decimation algorithms
* Boundary marking and other features
* Estimation of surface curvatures
* Generation of mesh surfaces around biological molecules

**Technical Features:**

* Code is implemented in C++ and supports Python using a pybind11 wrapper (pygamer).
* Cross system compilation using CMake and runs on Linux (64 bit), Windows (32 or 64 bit) or MacOS (64 bit).
* Blender addon which enables easy access to GAMer features using the pygamer API.
* Uses the Colored Abstract Simplicial Complex data ([CASC](http://github.com/ctlee/casc/) structure as the flexible underlying representation of surface and tetrahedral meshes.
* Code is hosted by [GitHub](http://github.com/ctlee/gamer/) under the Lesser GNU public license (LGPLv2). Please post issues or reports there.

## Acknowledging your use of GAMer
Thanks for using GAMer! The developers would love to hear how you are using the tool. Please send us an email or post on GitHub letting us know.

Please cite the above Zenodo DOI to acknowledge the software version and cite the following paper:<br/>
[Lee, C. T.; Laughlin, J. G.; Angliviel de La Beaumelle, N.; Amaro, R.; McCammon, J. A.; Ramamoorthi, R.; Holst, M. J.; Rangamani, P. GAMer 2: A System for 3D Mesh Processing of Cellular Electron Micrographs. bioRxiv 2019, 534479.](https://www.biorxiv.org/content/10.1101/534479v1)

## Installation
The following instructions are to build the base GAMer library.
If you wish to additionally compile the Blender GAMer addon, GAMer documentation, or other features please refer to the Additional Options section prior to building.

First, download a copy of the source from [releases](https://github.com/ctlee/gamer/releases) or clone the master branch.<br/>
```bash
git clone https://github.com/ctlee/gamer.git
cd gamer
```

Linux and Mac:
```bash
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr -DGAMER_TESTS=on -DCMAKE_BUILD_TYPE=Release ..
make
```

For Windows, we support building using Microsoft Visual Studio (MSVS) through the use of CMake generators:
```bash
mkdir build64
cd build64
cmake -DBUILD_BLENDGAMER=TRUE -G "Visual Studio 15 2017 Win64" -A x64 ..
cmake --build . --config Release
```

For a complete guide to installation, including configuration of PyGAMer and BlendGAMer please checkout the [online installation documentation](https://gamer.readthedocs.io/en/latest/install.html).

## External libraries bundled/downloaded with/by GAMer
* GAMer uses [Tetgen](http://wias-berlin.de/software/tetgen/) to generate
tetrahedralizations.

* GAMer uses [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) to
compute eigenvalues/eigenvectors of the Local Structure Tensor.

* GAMer uses [casc](https://github.com/ctlee/casc) as the underlying simplicial
complex data structure.

* GAMer uses [GoogleTest](https://github.com/google/googletest) to handle testing.

* GAMer uses [Pybind11](https://pybind11.readthedocs.io/en/stable/)

* Mesh checks in the GAMer Blender addon are inspired or borrowed from 3D Print Toolbox by Campbell Barton and Meshalyzer from CellBlender.

* [Triangle](https://www.cs.cmu.edu/~quake/triangle.html) is also bundled with GAMer but not currently used.


## Development Build Status

[![Development Build Status](https://travis-ci.org/ctlee/gamer.svg?branch=development)](https://travis-ci.org/ctlee/gamer)
[![Build status](https://ci.appveyor.com/api/projects/status/urffu7062fnohidl/branch/development?svg=true)](https://ci.appveyor.com/project/ctlee/gamer/branch/development)
[![Documentation Status](https://readthedocs.org/projects/gamer/badge/?version=development)](https://gamer.readthedocs.io/en/development?badge=development)
