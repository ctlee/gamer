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

* Code is implemented in C++ and supports Python using a ``pybind11`` wrapper (pygamer).
* Cross system compilation using ``CMake`` and runs on Linux (64 bit), Windows (32 or 64 bit) or MacOS (64 bit).
* Blender addon which enables easy access to ``GAMer`` features using the ``PyGAMer`` API.
* Uses the Colored Abstract Simplicial Complex data (`CASC <http://github.com/ctlee/casc/>`_) structure as the flexible underlying representation of surface and tetrahedral meshes.
* Code is hosted at `GitHub <http://github.com/ctlee/gamer/>`_ under the Lesser GNU public license (LGPLv2). Please post issues or reports there.


*******************************************
Acknowledging the Use of GAMer in Your Work
*******************************************

The `contributors <https://github.com/ctlee/gamer/contributors>`__ to this project are grateful for your use of this software.
To acknowledge our contributions please cite the following:

**Referencing the Code:** |GAMer Zenodo|

**Technical manuscript:**

.. code-block:: sh

  @article {Lee_gamer2_534479,
    author = {Lee, Christopher T. and Laughlin, Justin G.
              and Angliviel de La Beaumelle, Nils and Amaro, Rommie
              and McCammon, J. Andrew and Ramamoorthi, Ravi
              and Holst, Michael J. and Rangamani, Padmini},
    title = {GAMer 2: A System for 3D Mesh Processing of
             Cellular Electron Micrographs},
    elocation-id = {534479},
    year = {2019},
    doi = {10.1101/534479},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/early/2019/01/29/534479},
    eprint = {https://www.biorxiv.org/content/early/2019/01/29/534479.full.pdf},
    journal = {bioRxiv}
  }

.. |GAMer Zenodo| image:: https://zenodo.org/badge/122682242.svg
   :target: https://zenodo.org/badge/latestdoi/122682242


Legacy Contributions
====================

*   `GAMer 1 <http://fetk.org/codes/gamer/>`__ was originally developed by Zeyun Yu, Yuhui Cheng, and Michael Holst.
    To acknowledge your use of ``GAMer 1`` please cite:

    .. code-block::

      @article{Yu_gamer_2008,
        author = {Yu, Zeyun and Holst, Michael J.
                  and Cheng, Yuhui and McCammon, J. Andrew},
        title = {{Feature-Preserving Adaptive Mesh Generation for
                  Molecular Shape Modeling and Simulation}},
        journal = {J. Mol. Graph. Model.},
        issn = {10933263},
        year = {2008},
        month = jun,
        volume = {26},
        number = {8},
        pages = {1370--1380},
        doi = {10.1016/j.jmgm.2008.01.007},
        url = {https://dx.doi.org/10.1016/j.jmgm.2008.01.007},
      }

*   The ``BlendGAMer`` addon is inspired by work from Tom Bartol (Salk Institute) and Johan Hake who developed the original addon for ``GAMer 1``.

.. toctree::
   :maxdepth: 1
   :caption: Getting Started:

   install
   algorithms
   blendgamer
   faq
   changelog


.. toctree::
   :maxdepth: 2
   :caption: Tutorials and Examples:

   tutorial

.. toctree::
   :maxdepth: 1
   :caption: API Documentation:

   _cppapi/root
   pygamer

.. toctree::
   :maxdepth: 1
   :caption: For Developers:

   development

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
