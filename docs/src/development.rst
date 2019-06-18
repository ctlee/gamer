=================
Developer's Guide
=================

Contributing
============

The Basic Idea
--------------

We use the "development" branch to develop GAMer. When "development" has reached
a mature state in terms of current functionality and stability, we
merge "development" into the master branch. This happens at the time
when a release is made.

If you are developing or modifying features you should not work on "development"
directly. Rather, you should branch "development" into a new personal branch.
After testing your feature, you can submit a *pull request*. Your modifications
will be subject to testing by continuous integration and merged into the "development"
branch.


Issues, Bugs, and Feature Requests
----------------------------------

If you discover some unexpected behavior or wish to see a new feature in GAMer
please file an issue on `GitHub Issues`_.

.. _GitHub Issues: https://github.com/ctlee/gamer/issues

Members of the current development team will be able to advise or help you with
resolving the problem.

Testing
-------
Unit tests for the code are handled by pytest and google test for python and
C++ respectively. All tests can be run using ctest after a successful build:

.. code-block:: bash

    ctest -V

Documentation
-------------
Every function, class, and module that you write must be documented. Please check out

    https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt

how to do this. In short, after every function (class, module) header, there should
be a docstring enclosed by """ ... """, containing a short and long description what
the function does, a clear description of the input parameters, and the return value,
and any relevant cross-references or citations. You can include Latex-style math
in the docstring.

For a deeper understanding and reference please have a look at the Sphinx documentation

    http://sphinx-doc.org/

In particular, the API functions (those publicly visible to the external user) should
be well documented.

To build the documentation you need the dependencies from the file
``docs/requirements.txt`` which you can install via pip::

   pip install -r docs/rtd-requirements.txt

afterwards reconfigure your build with CMake you are ready to build the documentation::

   make sphinx_docs

The documentation can be found in ``{builddir}/docs/html/``.

Embedding Widgets in Jupyter-Notebooks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A notebook file may be saved with the current widget state as metadata. This allows the notebook file to be rendered with rendered widgets. To save a notebook with the current widget state, use the `Save Notebook Widget State` menu item.

In order to delete old saved state and save new state to the notebook, do the following in order:

1. Use the `Clear Notebook Widget State` menu and save the notebook. This clears the metadata from the notebook file.
2. Restart the kernel and refresh the page. This clears the old widget state from the widget manager on the page.
3. Create whatever widgets you'd like, and use `Save Notebook Widget State` and save the notebook. This saves the new widget state to the notebook file.


Release
-------

.. code-block:: bash
    [distutils]
    index-servers =
      pypi
      pypitest

    [pypi]
    username={VALUE}
    password={VALUE}

    [pypitest]
    repository=https://test.pypi.org/legacy/
    username={VALUE}
    password={VALUE}

`python setup.py sdist bdist_wheel`
`twine upload -r pypitest dist/*`
`pip install --index-url https://test.pypi.org/simple/ pygamer==0.0.14`
