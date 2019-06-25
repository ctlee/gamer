#########
Changelog
#########

******************
2.1.0 (??-??-????)
******************

**New Features**:

- Generalized decimation for Tetrahedral meshes (:pr:`33`). Thanks to :user:`ishantimalsina`.

******************
2.0.2 (??-??-????)
******************

**New Features**:

- Move from SWIG to use Pybind11 to generate Python binding (:pr:`30`).
- Now using FetchContent to populate external libraries (:pr:`30`).
- Introduced calculation of curvatures (:pr:`26`) and Helfrich energy evaluation (:pr:`27`).
- Substantially improved documentation of both C++ and Python hosted at https://gamer.readthedocs.io (:pr:`34`)

**Fixes**:

- Cleanup of compiler warnings and implicit casts (:pr:`30`).
- Various bug fixes (:pr:`19`, :pr:`24`, :pr:`29`). Thanks :user:`justinlaughlin`.

******************
2.0.1 (02-11-2019)
******************

**New Features**:

- Stable beta release! Compilation is supported on major operating systems. :pr:`16`