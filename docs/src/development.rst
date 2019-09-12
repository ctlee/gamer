*****************
Developer's Guide
*****************

************
Contributing
************

The Basic Idea
==============

We use the ``development`` branch to develop GAMer.
When ``development`` has reached a mature state in terms of current functionality and stability, we merge ``development`` into the ``master`` branch. This happens at the time when a release is made.

If you are developing or modifying features you should not work on ``development`` directly.
Rather, you should branch ``development`` into a new personal branch.
After testing your feature, you can submit a *pull request*. Your modifications will be subject to testing by continuous integration and merged into the ``development`` branch.


Issues, Bugs, and Feature Requests
==================================

If you discover some unexpected behavior or wish to see a new feature in GAMer please file an issue on `GitHub Issues`_.

.. _GitHub Issues: https://github.com/ctlee/gamer/issues

Members of the current development team will be able to advise or help you with resolving the problem.

Testing
=======

Unit tests for the code are handled by ``pytest`` and ``google test`` for Python and C++ respectively.
All tests can be run using ctest after a successful build:

.. code-block:: sh

    ctest -V

All new functionality introduced should be bundled with Unit Tests to ensure their successful operation and use.

Documentation
=============

Every function, class, and module that you write must be documented.
In short, every C++ function, namespace, and other should have a ``doxygen`` compatible description of what the function does, inputs, and outputs.
All Python functions should have a docstring added.
Note that github issues and users can be referenced using :issue:`XY` and :user:`username` markups.
This functionality is enabled by the `sphinx-issue` sphinx extension.

To build the documentation you need the dependencies from the file ``docs/requirements.txt`` which you can install via ``pip``.

.. code-block:: sh

   pip install -r docs/rtd-requirements.txt

afterwards reconfigure your build with CMake you are ready to build the documentation.

.. code-block:: sh

   make sphinx_docs

The documentation can be found in ``{builddir}/docs/html/``.


Embedding Widgets in Jupyter-Notebooks
--------------------------------------

A notebook file may be saved with the current widget state as metadata. This allows the notebook file to be rendered with rendered widgets.
To save a notebook with the current widget state, use the `Save Notebook Widget State` menu item.

In order to delete old saved state and save new state to the notebook, do the following in order:

#. Use the `Clear Notebook Widget State` menu and save the notebook. This clears the metadata from the notebook file.
#. Restart the kernel and refresh the page. This clears the old widget state from the widget manager on the page.
#. Create whatever widgets you'd like, and use `Save Notebook Widget State` and save the notebook. This saves the new widget state to the notebook file.

Inline Code Comments
--------------------

You are also encouraged to use traditional comments inline with the code to explain complex steps.
We also encourage the use of ``TODO`` style comments in the source.
These tags can be parsed by tools such as the `TodoReview <https://packagecontrol.io/packages/TodoReview>`__ plugin for Sublime Text to generate useful roadmaps.
Comment tags should follow this style:

.. code-block:: cpp

   // TODO: (#priority) {Text comments}
   // TODO: (100) Example todo note

Where ``#priority`` is a numerical priority ranging from 0 to 100 for sorting importance with 0 being the most important.


Publishing a new official release
=================================

#.  merge the `development` branch into `master`.

    .. code-block:: sh

      git checkout master; git merge development

#.  make a new tag 'v{major}.{minor}.{patch}'. To determine the new version follow the guidelines outlined by `Semantic Versioning <https://semver.org/>`__.

    .. code-block:: sh

      git tag -a v2.0.1 -m "Description of the release"

#.  Push the new tag to the remote repository.

    .. code-block:: sh

      git push origin v2.0.1

#.  Update PyPi distribution. First test the distribution package accordingly.

    .. code-block:: bash

      python setup.py sdist bdist_wheel
      twine upload -r pypitest dist/*
      pip install --index-url https://test.pypi.org/simple/ pygamer==0.0.14

    It may be helpful to declare `export PIP_NO_BUILD_ISOLATION=false` since many projects are not available on the test PyPi server.
    In accord with PEP518 and PEP517, pip will attempt to grab build depdencies in isolation and will throw errors when a required library cannot be found.