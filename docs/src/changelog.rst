#########
Changelog
#########

******************
2.1.0 (??-??-????)
******************

**New Features**:

- Generalized decimation for Tetrahedral meshes (:pr:`33`). Thanks to :user:`ishantimalsina` and :user:`mvhsan`.

******************
2.0.4 (09-12-2019)
******************

**New Features**:

- Updated `BlendGAMer` to work with newly released `Blender` 2.80.
- Improvements to the documentation.

**Fixes**:

- Fixes for matplotlib dependency detection in Blender Addon.

******************
2.0.3 (07-24-2019)
******************

**New Features**:

- Faster mesh curvature calculation with better Blender interface (:pr:`38`).
- Edge flip now considered the vertex adjacency (:pr:`38`).
  We now impose a penalty/boost to prefer reasonable adjacencies.
- Normal vectors of faces and vertices can now be cached to improve speed (:pr:`38`).

**Fixes**:

- Fix for Dolfin XML writer. Previously vertex renumbering in the writer could generate mislabeled boundary markings.
  It is recommended that all users of this function update to this release.

******************
2.0.2 (06-25-2019)
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

- Stable beta release! Compilation is supported on major operating systems (:pr:`16`).