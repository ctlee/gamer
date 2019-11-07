##########################
Frequently Asked Questions
##########################

#.  How do I resolve the ``ModuleNotFoundError: No module named 'skbuild'`` when using ``pip`` to install ``PyGAMer``?

    The name of the package is ``scikit-build`` which unfortunately does not match the module name.
    Run ``pip install scikit-build`` to satisfy the requirement.