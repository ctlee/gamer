###########################################
GAMer - Geometry-preserving Adaptive Mesher
###########################################

GAMer is a surface mesh improvement library developed to condition surface meshes derived from noisy biological imaging data.
Using Tetgen, GAMer can
generate tetrahedral meshes suitable for finite elements simulations of
reaction-diffusion systems among other applications.
GAMer has the following main features:

* Surface mesh improvement and decimation algorithms
* Boundary marking and other features
* Estimation of surface curvatures
* Generation of mesh surfaces around biological molecules

**Technical Features:**

* Code is implemented in C++ and supports Python using a pybind11 wrapper (pygamer).
* Cross system compilation using CMake and runs on Linux (64 bit), Windows (32 or 64 bit) or MacOS (64 bit).
* Blender addon which enables easy access to GAMer features using the pygamer API.
* Uses the Colored Abstract Simplicial Complex data (`CASC <http://github.com/ctlee/casc/>`_) structure as the flexible underlying representation of surface and tetrahedral meshes.
* Code is hosted at `GitHub <http://github.com/ctlee/gamer/>`_ under the Lesser GNU public license (LGPLv2). Please post issues or reports there.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started:

   install
   changelog

.. toctree::
   :maxdepth: 2
   :caption: Tutorials and Examples:

   tutorial

.. toctree::
   :maxdepth: 2
   :caption: API Documentation:

   _cppapi/root
   pygamer

.. toctree::
   :maxdepth: 2
   :caption: For Developers:

   development

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
