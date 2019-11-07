############
Installation
############

GAMer is first and foremost a C++ library.
If you are an experienced C++ developer and wish to link to ``libGAMer`` start at :ref:`Building libGAMer`.
For other users, many of the functions and types provided by GAMer can be accessed through the Python API wrapper (``PyGAMer``) or can be called through an interactive Blender plugin (``BlendGAMer``).
Check out the :ref:`Getting PyGAMer` and :ref:`Getting BlendGAMer` sections for details on how to install these tools on your system.

.. contents:: Skip to a Section
   :local:

.. _Prerequisites:

*************
Prerequisites
*************

In order to build GAMer, ``PyGAMer``, and ``BlendGAMer``, you will need several tools and libraries on your system.

  * ``C++ compiler``: supporting C++14 standard or newer.
  * ``CMake``: version 3.10 or newer.
  * ``Eigen 3``: CMake can download and locally configure Eigen for you.
  * ``CASC``: CMake can download and locally configure CASC for you.
  * ``pybind11``: (Required for ``PyGAMer`` and ``BlendGAMer``) CMake can download and locally configure PyBind11 for you.
  * ``Python Headers``: (Required for ``PyGAMer`` and ``BlendGAMer``) This is often available by default on your system. Alternatively Python environments can be configured using `Anaconda <https://www.anaconda.com/>`__.

  * ``scikit-build``: Required for build and install using ``pip`` and ``setup.py``. Obtain using ``pip install scikit-build``.

  * ``Blender``: (Optional) To build the Blender addon, it is preferable to have a working installation of `Blender <https://www.blender.org/>`__ on your path. This way CMake can verify addon compatibility prior to packaging.

  * ``Doxygen``: (Optional) is necessary to build the standalone C++ API documentation.
  * ``Breathe``, ``Exhale``, ``Sphinx``: (Optional) Python libraries are used to build the Python/C++ API documentation. Other Python dependencies such as those necessary for building are handled by ``pip`` and ``setuptools``.

.. _Building libGAMer:

*****************
Building libGAMer
*****************

The provided steps will compile **both the shared and static libraries**.
For users of CMake 3.12 or later, there is very little overhead to compiling both libraries since we are taking advantage of CMake's object library capabilities to compile the sources only once.

If you wish to additionally compile the Blender GAMer addon, GAMer documentation, or other features please refer to the :ref:`Additional CMake Options` section prior to building.
If you are compiling on Windows using Microsoft Visual Studio please also refer to the :ref:`Compiling on Windows` section.

#.  Download a copy of the source from `GitHub Releases <https://github.com/ctlee/gamer/releases>`__ or by cloning the master/development branch.

    .. code-block:: sh

      git clone https://github.com/ctlee/gamer.git

#.  For a traditional build and installation, create an out of source folder.

    .. code-block:: sh

      cd gamer
      mkdir build; cd build

#.  Run CMake to configure and build the project.

    .. code-block:: sh

      cmake -DGAMER_TESTS=on -DCMAKE_BUILD_TYPE=Release ..
      cmake --build . --config Release -j 2

#.  Run the unit tests to ensure successful compilation.

    .. code-block:: sh

      ctest -C Release -V

#.  Install the libraries and headers.

    .. code-block:: sh

      cmake --build . --target install

.. _Alternative Python Build:

Alternative Python Build
========================

Alternatively you can use ``setup.py`` which is configured to use ``scikit-build`` and interfaces with CMake to build the library.
This should work without modification on all supported platforms.

.. code-block:: sh

  python setup.py build
  python setup.py install

This performs approximately the same CMake build steps as above, however in an automated fashion.
You can pass :ref:`Additional CMake Options` through the setup by appending them to the setup call.
Other details about the scikit-build process can be found `here <https://scikit-build.readthedocs.io/en/latest/>`__.

.. code-block:: sh

  python setup.py build -- -DCMAKEOPT=...

.. Warning::

  scikit-build will install package components to Python specific locations. You may need to adjust your includes search paths to help your compiler find the relevant GAMer header and library files.


.. _Compiling on Windows:

Compiling on Windows
====================

For Windows, we support building using Microsoft Visual Studio (MSVS).
The process is essentially the same as the traditional build except that the `CMake MSVS generator`_ expects an architecture.
To build a 64-bit library you need only append ``-A x64`` to the initial CMake configure.
You can also use the alternative Python build process with no additional modifications as it detects your Python bit version and matches it.

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
  * -  Build the ``PyGAMer`` extension.
    - ``-DBUILD_PYGAMER=on``
  * - Specify the Python executable path.
    - ``-DPYTHON_EXECUTABLE:FILEPATH=/path/to/python3``
  * - Package the Blender addon. This flag automatically builds the Python extension.
    - ``-DBUILD_BLENDGAMER=on``
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
  * - Install ``BlendGAMer`` to the user Blender addon path. This requires Blender to be on your systems PATH.
    - ``-DBLENDER_PLUGIN_INSTALL=on``
  * - Enforce strict Python version matching with Blender.
    - ``-DBLENDER_VERSION_STRICT=on``
  * - Compile the Tetgen binary.
    - ``-DBUILD_TETGEN_BIN=on``
  * - Specify the installation prefix for GAMer headers and libraries
    - ``-DCMAKE_INSTALL_PREFIX=/usr/local``

.. _Getting PyGAMer:

***************
Getting PyGAMer
***************

.. _PyGAMer the Easy Way:

The Easy Way
============

.. note::

   ``PyGAMer`` is developed for use with Python 3 and newer.
   Other Python versions may work but are untested and may require workarounds.

We recommend that you install ``PyGAMer`` using the ``pip``>=10.0 utility.

.. code-block:: sh

  pip install pygamer

The pip utility will automatically sort out the package dependencies for you and potentially build the library.
Unfortunately ``pip`` is not traditionally bundled with the prepackaged Blender installation consult the :ref:`Getting BlendGAMer` section for instructions on how to build ``BlendGAMer``.

.. warning::

   ``pip`` versions < 10.0 do not implement PEP 518 which specifies the installation and local containment of build time dependencies.
   These can be be resolved by installing these dependencies from a `requirements.txt <https://github.com/ctlee/gamer/blob/master/requirements.txt>`__ file prior to running ``pip install pygamer``.

.. _PyGAMer the Harder Way:

The Harder Way
==============

.. note::

   ``PyGAMer`` is developed for use with Python 3 and newer.
   Other Python versions may work but are untested and may require workarounds.

You can also build ``PyGAMer`` using ``setuptools`` on your own using the alternative build instructions.
By default, the Python setup enables the compilation of the ``PyGAMer`` Python extension module.

.. code-block:: sh

  python setup.py build
  python setup.py install

If you insist on it, it is also possible to build and install ``PyGAMer`` using CMake directly.
This will place the plugin into ``${PYTHON_SITE_PACKAGES}/pygamer/*``.
Although the CMake Python module installation is available, it can be error prone and therefore we recommend building using ``setuptools``.

.. code-block:: sh

  mkdir build; cd build
  cmake -DBUILD_PYTHONEXT=on ..
  cmake --build . --config Release -j 2
  cmake --build . --target install


.. _Getting BlendGAMer:

******************
Getting BlendGAMer
******************

.. _BlendGAMer the Easy Way:

.. warning::
   Currently ``BlendGAMer`` only supports ``Blender`` v2.79b.
   If you have another version of ``Blender``, please install `Blender v2.79b <https://download.blender.org/release/Blender2.79/>`__ before proceeding.

The Easy Way
============

If you seek to use GAMer in an interactive way through Blender.
The easiest way to get a working copy is to download and install one of the prebuilt binaries of the GAMer Blender addon are available under `github releases`_.
These zip files contain prebuilt ``PyGAMer`` binaries which correspond to specific Blender release versions published by the Blender Foundation.
The zip file can be installed by following the traditional `Blender addon installation instructions`_.

.. Warning::
  If you are using a non-standard installation, such as Blender you have compiled yourself or from a package distribution (i.e., ``apt`` or ``yum``), the precompiled zip addons may not work for you.
  This is because the Python extension module version must be compiled using Python version matching Blender's bundled Python version.
  Package distributions often use Python versions already available on your system and therefore the precompiled binaries may not match.
  To resolve this, you will need to install :ref:`BlendGAMer the Harder Way`.

.. _github releases: https://github.com/ctlee/gamer/releases

.. _Blender addon installation instructions: https://docs.blender.org/manual/en/latest/preferences/addons.html#header


.. _BlendGAMer the Harder Way:

The Harder Way
==============

You can build ``BlendGAMer`` yourself using ``CMake``.
Owing to the complexities of building Python extension modules, it is preferable to have a working installation of Blender on your system.
While this is not strictly necessary, it enables CMake to verify that the Python versions will be compatible.
Note that the prebuilt Blender binaries from the Blender Foundation do not contain Python header files and are therefore unsuitable for compilation.

#.  Ensure you have a working Blender installation. And if possible append the Blender executable to your systems ``PATH``.
    Follow instructions online for `Getting Blender <https://docs.blender.org/manual/en/latest/getting_started/installing/>`__.

    On Mac add the following commands to ``~/.bash_profile`` pointing to the directory with Blender's binary:

    .. code-block:: sh

      export PATH="/Applications/blender/blender.app/Contents/MacOS:${PATH}"

    On Linux add the following command to ``~/.bashrc`` or ``~/.profile`` pointing to the directory with Blender's binary:

    .. code-block:: sh

      export PATH=/path/to/blender/directory/bin:$PATH

    On Windows, execute from the command line:

    .. code-block:: sh

      blender -r

#.  Check if your version of Blender is bundled with its own Python.

    .. code-block:: sh

      blender -b --factory-startup --python-expr "import bpy; print(bpy.app.binary_path_python);"

    - A)  If the printed string indicates a Python binary inside of a Blender folder e.g., ``/Applications/blender/blender.app/Contents/Resources/2.79/python/bin/python3.5m`` or ``/usr/local/blender/2.79b/2.79/python/bin/python3.5m``, this indicative of a bundled Python.

    - B)  If the Python binary path is not under a Blender folder e.g., ``/usr/bin/python3.6m`` then Blender is using some other Python distribution.

#.  Setup your Python development environment.

    - A)  For Blender with bundled Python you will need to get a separate Python development environment.
          To get a compatible Python suitable for building ``BlendGAMer`` we recommend using `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__.
          Check the table below to verify the Python version for your target Blender version.

          ===============  ==============
          Blender Version  Python Version
          ===============  ==============
          2.79b            3.5
          ===============  ==============

          Create a new environment corresponding to the Python version.

          .. code-block:: sh

             conda create --name py35 python=3.5
             conda activate py35

          .. note:: Don't forget to activate your new Python environment if you open a new shell.

    - B)  Depending on your configuration you may need to install the development headers and ``numpy``. For example:

          .. code-block:: sh

             sudo apt install python3.6-dev python3-numpy python3-pip

#.  Configure, and build ``libGAMer``, ``PyGAMer``, and ``BlendGAMer``.

    .. code-block:: sh

       cmake -DBUILD_BLENDGAMER=on -DGAMER_TESTS=on -DCMAKE_BUILD_TYPE=Release ..
       cmake --build . --config Release -j 2


#.  Install! At this point you should have a packaged ``.zip`` at the root of your out-of-source build directory.
    Follow the `Blender addon installation instructions`_ to install.

    Alternatively you can have CMake install the addon into Blender's User addons folder.

    .. code-block:: sh

       cmake -DBLENDER_PLUGIN_INSTALL=on ..
       cmake --build . --target install

#.  Load up Blender and verify that ``BlendGAMer`` is working maybe by following one of our illustrative :ref:`BlendGAMer Tutorials`.

.. _Getting matplotlib in Blender:

Installing Matplotlib in Blender
================================

For advanced users only, if you wish to run curvature calculations in ``BlendGAMer`` there is a ``matplotlib`` dependency which is not satisfied by default ``Blender``.
For ``Blender`` versions using the bundled system ``Python``, you may only need to install the relevant ``python3-matplotlib`` or related package for your system.

Otherwise, if you are using a prepackaged version of ``Blender``, the currently recommended method to get ``matplotlib`` is through ``pip``.
Fist download the ``get-pip.py`` file from the `pip documentation <https://pip.pypa.io/en/stable/installing/>`__.
Execute this script using the bundled ``Python`` from ``Blender``.

.. code-block:: sh

   python get-pip.py

The bundled ``Python`` can be found at

.. code-block:: sh

   {path to blender}/2.xx/python/bin/python

for Linux and Windows and at

.. code-block:: sh

   /blender.app/Contents/Resources/2.79/python/bin

for Mac platforms.
Now that ``pip`` is installed you can use it to install ``matplotlib``:

.. code-block:: sh

   /path/to/blenderspython/pip install matplotlib

``matplotlib`` should now be installed.


.. _Building the Documentation:

**************************
Building the Documentation
**************************

You can always read the latest documentation online on `Read The Docs <https://gamer.readthedocs.io>`__.
The preferred way to build the documentation is through the use of ``setuptools``.
Be sure to append ``-DGAMER_DOCS=on`` to your call to download the external repository of tutorials.

.. code-block:: sh

    python setup.py install -- -DGAMER_DOCS=on
    python setup.py build_sphinx

Other options to ``setuptools`` can be found at `Sphinx setuptools integration <https://www.sphinx-doc.org/en/master/usage/advanced/setuptools.html>`__.

It is also possible to compile the documentation using CMake by building target ``sphinx_docs``.

.. note::
    If you are getting a module import error, this is indicative that Python cannot find an installed copy of ``PyGAMer`` to retrieve docstrings from.
    You can either manually append the location of the ``PyGAMer`` extension module to the PYTHONPATH in ``docs/conf.py.in``.
    Alternatively you can install ``PyGAMer`` in a more conventional location.