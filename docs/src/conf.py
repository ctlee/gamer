# Config to simplify building on ReadTheDocs!
# The Sphinx documentation can also be built locally using CMake.
#
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../../build/lib/'))

# -- Project information -----------------------------------------------------

project = 'GAMer'
copyright = '2019, Christopher T. Lee'
author = 'Christopher T. Lee'

# The full version, including alpha/beta/rc tags
# release = '2.0.1 beta'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
'sphinx.ext.autosummary',
'sphinx.ext.napoleon',
'sphinx.ext.intersphinx',
'breathe',
'exhale'
]

autosummary_generate = True

breathe_projects = { "gamer_project": "./_doxyoutput/xml" }
breathe_default_project = "gamer_project"

# Setup the exhale extension
exhale_args = {
    ############################################################################
    # These arguments are required.                                            #
    ############################################################################
    "containmentFolder":     "./_cppapi",
    "rootFileName":          "root.rst",
    "rootFileTitle":         "C++ API Reference",
    "doxygenStripFromPath":  "..",
    ############################################################################
    # Suggested optional arguments.                                            #
    ############################################################################
    "createTreeView":        True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin":    "INPUT = ../../include",
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '_templates']


# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# html_static_path = ['/Users/ctlee/gamer/docs/_static']

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'https://docs.python.org/': None}
