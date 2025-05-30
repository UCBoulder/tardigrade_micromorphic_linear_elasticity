.. _changelog:


#########
Changelog
#########

******************
1.6.0 (unreleased)
******************

******************
1.5.1 (2025-05-21)
******************

New Features
============
- Allow the user to perform a full stack build of tardigrade (:pull:`23`). By `Nathan Miller`_.

Internal Changes
================
- Preparing for conda packaging (:pull:`24`). By `Nathan Miller`_.

Release
=======
- Released version (:pull:`25`). By `Nathan Miller`_.

******************
1.5.0 (2024-11-26)
******************

Release
=======
- Released version (:pull:`22`). By `Nathan Miller`_.

Breaking Changes
================
- Removed the old subroutines since hydra is the preferred method of use (:pull:`20`). By `Nathan Miller`_.

Documentation
=============
- Updated changelog to reflect changes to version (:pull:`21`). By `Nathan Miller`_.

******************
1.4.0 (2024-11-25)
******************

Release
=======
- Released version (:pull:`19`). By `Nathan Miller`_.

New Features
============
- Added additional information to output messages in case of errors (:pull:`13`). By `Nathan Miller`_.

Internal Changes
================
- Removed comma causing warnings in cmake builds (:pull:`7`). By `Nathan Miller`_.
- Use new error_tools format and allow for timestep cutbacks (:pull:`10`). By `Nathan Miller`_.
- Added flag that will allow all errors to be treated as convergence errors (:pull:`11`). By `Nathan Miller`_.
- Enabled use of tardigrade_hydra's preconditioner (:pull:`12`). By `Nathan Miller`_.
- Allowing for changes to hydra base class inputs (:pull:`14`). By `Nathan Miller`_.
- Allow for the package version to be specified when using FetchContent (:pull:`16`). By `Nathan Miller`_.

Bug Fixes
=========
- Enabled compilation with the new tardigrade_micromorphic_tools (:pull:`8`). By `Nathan Miller`_.
- Removed trailing whitespace from add_library (:pull:`15`). By `Nathan Miller`_.

******************
1.3.2 (2024-01-24)
******************

Release
=======
- Released version (:pull:`6`). By `Nathan Miller`_.

New Features
============
- Added a hydra-based micromorphic linear elasticity (:pull:`3`). By `Nathan Miller`_.

Internal Changes
================
- Enabled github actions (:pull:`1`). By `Nathan Miller`_.
- Added tardigrade_hydra as a build requirement (:pull:`2`). By `Nathan Miller`_.
- Minor update to documentation (:pull:`4`). By `Nathan Miller`_.
- Updated changelog prior to release (:pull:`5`). By `Nathan Miller`_.

******************
1.3.1 (2023-07-25)
******************

Breaking changes
================
- Change project, package, and namespace to use the 'tardigrade' prefix (:issue:`7`, :merge:`18`). By `Kyle Brindley`_.

******************
1.2.1 (2023-07-11)
******************

Internal Changes
================
- Replace shell scripts with inline Gitlab CI configuration (:issue:`2`, :merge:`12`). By `Kyle Brindley`_.
- Create project specific CI environment (:issue:`3`, :merge:`13`). By `Kyle Brindley`_.
- Use setuptools_scm for Git tag versioning (:issue:`4`, :merge:`14`). By `Kyle Brindley`_.
- Conda package and deployment (:issue:`5`, :merge:`15`). By `Kyle Brindley`_.

Bug Fix
=======
- Identified and corrected access out of bounds error in test of gradient of the reference stresses w.r.t.
  the gradient of the micro deformation (:merge:`17`). By `Nathan Miller`_.

******************
1.1.1 (2022-11-03)
******************

Internal Changes
================

- Fixed bug in linear elastic constraint equations (:merge:`7`). By `Nathan Miller`_.
- Fixed additional bug in linear elastic constraint equations (:merge:`8`). By `Nathan Miller`_.

******************
1.1.0 (08-16-2022)
******************

Internal Changes
================

- Moved the code to the cpp_stub format (:merge:`1`). By `Nathan Miller`_.
- Moved the tests to the BOOST test format (:merge:`2`). By `Nathan Miller`_.
- Removed old material library interface definitions (:merge:`3`). By `Nathan Miller`_.
- Added the ability to turn of building the python bindings (:merge:`4`). By `Nathan Miller`_.
- Added wrapper for calculation of current stresses from the fundamental deformation measures (:merge:`5`). By `Nathan Miller`_.

Release
=======

- Released version 1.1.0 (:merge:`6`). By `Nathan Miller`_.
