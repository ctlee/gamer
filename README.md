# Release snapshot of GAMer
GAMer is a surface mesh improvement library included in the FEtk 
software umbrella.

## Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites
To use this library you will need to a C++ compiler which has full support of C++14.
Note that there is a known issue with Clang 4.x.x versioned compilers (including XCode version 9), which will produce a compiler error. 
The current workaround to this problem is to obtain and use Clang version 5+.

### Installing
```bash
git clone https://github.com/ctlee/gamer.git
cd gamer
git submodule init
git submodule update
```

```bash
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr -DBUILD_TESTS=on -DCMAKE_BUILD_TYPE=Release ..
make -j12
```

## Authors
* **John Moody** <jbmoody@ucsd.edu>
* **Christopher Lee** <ctlee@ucsd.edu>

## License
This project is licensed under the GNU Lesser General Public License v2.1 - 
please see the [COPYING.md](COPYING.md) file for details.

## Acknowledgments
This project is supported by the National Institutes of Health under grant
numbers P41-GM103426, T32-GM008326, and R01-GM31749. 
It is also supported in part by the National Science Foundation under awards DMS-CM1620366 and DMS-FRG1262982.


