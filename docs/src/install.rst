############
Installation
############

GAMer is first and foremost a C++ library.
If you are an experienced C++ developer and wish to link to libGAMer start at :ref:`Building libGAMer`.
For other users, many of the functions and types provided by GAMer can be accessed through the Python API wrapper (PyGAMer) or can be called through an interactive Blender plugin (BlendGAMer).
Check out the :ref:`Getting PyGAMer` and :ref:`Getting BlendGAMer` sections for details on how to install these tools on your system.


.. _Prerequisites:

*************
Prerequisites
*************

In order to build GAMer, PyGAMer, and BlendGAMer, you will need several tools and libraries on your system.
Other than the build toolchain and

  * ``C++ compiler``: supporting C++14 standard or newer.
  * ``CMake``: version 3.11 or newer.
  * ``Eigen 3``: CMake can download and locally configure Eigen for you.
  * ``CASC``: CMake can download and locally configure CASC for you.
  * ``pybind11``: (Required for PyGAMer and BlendGAMer) CMake can download and locally configure PyBind11 for you.
  * ``Python Headers``: (Required for PyGAMer and BlendGAMer) This is often available by default on your system. Alternatively Python environments can be configured using `Anaconda <https://www.anaconda.com/>`__.

  * ``Blender``: (Optional) To build the Blender addon, it is preferable to have a working installation of `Blender <https://www.blender.org/>`__ on your path. This way CMake can verify addon compatibility prior to packaging.

  * ``Doxygen``: (Optional) is necessary to build the standalone C++ API documentation.
  * ``Breathe``, ``Exhale``, ``Sphinx``: (Optional) Python libraries are used to build the Python/C++ API documentation. Other Python dependencies such as those necessary for building are handled by pip and setuptools.

.. _Building libGAMer:

*****************
Building libGAMer
*****************

The provided steps will compile **both the shared and static libraries** by default.
There is very little overhead to compiling both libraries since we are taking advantage of CMake's object library capabilities to compile the sources only once.

If you wish to additionally compile the Blender GAMer addon, GAMer documentation, or other features please refer to the :ref:`Additional CMake Options` section prior to building.
If you are compiling on Windows using Microsoft Visual Studio please also refer to the :ref:`Compiling on Windows` section.

1. Download a copy of the source from `GitHub Releases <https://github.com/ctlee/gamer/releases>`__ or by cloning the master/development branch.

::

  git clone https://github.com/ctlee/gamer.git

2. For a traditional build and installation, create an out of source folder.

::

  cd gamer
  mkdir build; cd build

3. Run CMake to configure and build the project.

::

  cmake -DGAMER_TESTS=on -DCMAKE_BUILD_TYPE=Release ..
  cmake --build . --config Release -j 2

4. Run the unit tests to ensure successful compilation.

::

  ctest -C Release -V

5. Install the libraries and headers.

::

  cmake --build . --target install


Alternative Python build
========================

Alternatively you can use ``setup.py`` which is configured to use ``scikit-build`` and interfaces with CMake to build the library.
This should work without modification on all supported platforms.

::

  python setup.py build
  python setup.py install

This performs approximately the same CMake build steps as above, however in an automated fashion.
You can pass :ref:`Additional CMake Options` through the setup by appending them to the setup call.
Other details about the scikit-build process can be found `here <https://scikit-build.readthedocs.io/en/latest/>`__.

::

  python setup.py build -- -DCMAKEOPT=...

.. Warning::

  scikit-build will install package components to Python specific locations. You may need to adjust your includes search paths to help your compiler find the relevant GAMer header and library files.


.. _Compiling on Windows:

Compiling on Windows
====================

For Windows, we support building using Microsoft Visual Studio (MSVS).
The process is essentially the same as the traditional build except that the `CMake MSVS generator`_ expects an architecture.
To build a 64-bit library you need only append ``-A x64`` to the initial CMake configure.
You can also use the alternative Python build process with no additional modifications.

.. _CMake MSVS generator: https://cmake.org/cmake/help/latest/generator/Visual%20Studio%2015%202017.html

::

  mkdir build64
  cd build64
  cmake -DGAMER_TESTS=on -A x64 ..
  cmake --build . --config Release -j 2

.. note::

  If you get an "ImportError: DLL load failed" you are likely linking a
  different python library version than Blender's bundled python.
  We recommend using Anaconda to obtain a python version matching Blender.


.. _Additional CMake Options:

Additional CMake Options
========================

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
    - ``-DSINGLE=on``
  * - Download the external GAMer documentation.
    - ``-DGAMER_DOCS=on``
  * - Configure the test cases.
    - ``-DGAMER_TESTS=on``
  * - Verbose configuration.
    - ``-DGAMER_CMAKE_VERBOSE=on``
  * - Download pybind11 locally
    - ``-DGETPYBIND11=on``
  * - Download Eigen 3 locally
    - ``-DGETEIGEN=on``

**Special options:**

.. list-table::
  :widths: 50 50
  :header-rows: 1

  * - Explanation
    - CMake Directive
  * - Install BlendGAMer to the user Blender addon path. This requires Blender to be on your systems PATH.
    - ``-DBLENDER_PLUGIN_INSTALL=on``
  * - Enforce strict Python version matching with Blender.
    - ``-DBLENDER_VERSION_STRICT=on``
  * - Compile the Tetgen binary.
    - ``-DBUILD_TETGEN_BIN=on``

.. _Getting PyGAMer:

***************
Getting PyGAMer
***************

The Easy Way
============

We recommend that you install PyGAMer using the pip utility.

::

  pip install pygamer

The pip utility will automatically sort out the package dependencies for you and potentially build the library.
Unfortunately pip is not traditionally bundled with the prepackaged Blender installation consult the :ref:`Getting BlendGAMer` section for instructions on how to build BlendGAMer.

The Harder Way
==============

You can also build PyGAMer using setuptools on your own using the alternative build instructions.
By default, the Python setup enables the compilation of the PyGAMer Python extension module.

::

  python setup.py build
  python setup.py install

If you insist on it, it is also possible to build and install PyGAMer using CMake directly.
This will place the plugin into ``${PYTHON_SITE_PACKAGES}/pygamer/*``.
Although the CMake Python module installation is available, it can be error prone and therefore we recommend building using setuptools and scikit-build.

::

  mkdir build; cd build
  cmake -DBUILD_PYTHONEXT=on ..
  cmake --build . --config Release -j 2
  cmake --build . --target install


.. _Getting BlendGAMer:

******************
Getting BlendGAMer
******************

The Easy Way
============

If you seek to use GAMer in an interactive way through Blender. The easiest way
to get a working copy is to download and install one of the prebuilt binaries
of the GAMer Blender addon are available under `github releases`_. The zip file
can be installed by following the traditional `Blender addon installation instructions`_.

.. _github releases: https://github.com/ctlee/gamer/releases

.. _Blender addon installation instructions: https://docs.blender.org/manual/en/latest/preferences/addons.html#header

The Harder Way
==============

You can build BlendGAMer yourself using CMake to help...


.. _Building the Documentation:

**************************
Building the Documentation
**************************

The preferred way to build the documentation is through the use of setuptools.

.. code-block:: bash

    python setup.py install -- -DGAMER_DOCS=on
    python setup.py build_sphinx

It is also possible to compile the documentation using CMake by building target ``sphinx_docs``.

.. note::
    If you are getting a module import error, this is indicative that Python cannot find an installed copy of PyGAMer to retrieve docstring from.
    You can either manually append the location of the PyGAMer extension module to the PYTHONPATH in ``docs/conf.py.in``.
    Alternatively you can install PyGAMer in a more conventional location.