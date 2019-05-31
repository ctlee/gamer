**********************
Prebuilt Blender Addon
**********************

If you seek to use GAMer in an interactive way through Blender. The easiest way
to get a working copy is to download and install one of the prebuilt binaries
of the GAMer Blender addon are available under `github releases`_. The zip file
can be installed by following the traditional `Blender addon installation instructions`_.

.. _github releases: https://github.com/ctlee/gamer/releases

.. _Blender addon installation instructions: https://docs.blender.org/manual/en/latest/preferences/addons.html#header

*********************
Compiling from Source
*********************

Prerequisites
-------------

In order to build GAMer you will need the following libraries available on your
system:

* C++ compiler supporting C++14 standard or newer.
* CMake version 3.11 or newer.
* (Optional) Eigen 3 library. If a copy cannot be found, CMake will automatically
  download it for you.
* (Optional) To build the Blender addon, it is prefereable to have a working
  installation of Blender. This way CMake can verify library compatibility prior
  to packaging.

.. _Configuring and Building:

Configuring and Building
------------------------

The following instructions are to build the base GAMer library.
If you wish to additionally compile the Blender GAMer addon, GAMer
documentation, or other features please refer to the :ref:`Additional CMake Options`
section prior to building. If you are compiling on Windows using Microsoft
Visual Studio please also refer to the :ref:`Compiling on Windows` section.

First, download a copy of the source by cloning the master branch.

.. code-block:: bash

    git clone https://github.com/ctlee/gamer.git
    cd gamer
    mkdir build
    cd build
    cmake -DGAMER_TESTS=on -DCMAKE_BUILD_TYPE=Release ..
    make

.. _Compiling on Windows:

Compiling on Windows
--------------------

For Windows, we support building using Microsoft Visual Studio (MSVS) through the use of CMake generators:

.. code-block:: bash

    mkdir build64
    cd build64
    cmake -DBUILD_BLENDER=TRUE -A x64 ..
    cmake --build . --config Release -j 2

.. note::
    If you get an "ImportError: DLL load failed" you are likely linking a
    different python library version than Blender's bundled python.
    We recommend using Anaconda to obtain a python version matching Blender.


.. _Additional CMake Options:

Additional CMake Options
------------------------
To enable these additional options append the flags to your initial CMake function call.
These can be used in addition to the standard `CMake flags`_.

.. _CMake flags: https://cmake.org/cmake/help/latest/manual/cmake.1.html

.. list-table::
    :widths: 50 50
    :header-rows: 1

    * - Explanation
      - CMake Directive
    * -  Build the pygamer extension.
      - ``-DBUILD_PYGAMER=on``
    * - Specify the Python executable path.
      - ``-DPYTHON_EXECUTABLE:FILEPATH=/path/to/python3``
    * - Package the Blender addon. This flag automatically builds the Python extension.
      - ``-DBUILD_BLENDER=on``
    * - Use single precision floating point numbers.
      - ``-DSINGLE=ON``
    * - Compile the Tetgen binary.
      - ``-DBUILD_TETGEN_BIN=on``
    * - Configure the test cases.
      - ``-DGAMER_TEST=on``
