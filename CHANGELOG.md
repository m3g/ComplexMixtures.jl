ComplexMixtures.jl Changelog
===========================
  
[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-experimental]: https://img.shields.io/badge/Experimental-yellow.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-fix]: https://img.shields.io/badge/Fix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg

Version 2.14.4
-------------
- ![ENHANCEMENT][badge-enhancement] Read chemfiles trajectories in place.

Version 2.14.3
-------------
- ![BUGFIX][badge-bugfix] fix diagonal unitcell test for when cell is not orthorhombic and contains negative vector entries.
- ![INFO][badge-info] Update FortranFiles dependency to 0.6.2 and use `seekstart` instead of `rewind`. 

Version 2.14.2
-------------
- ![ENHANCEMENT][badge-enhancement] Better finalizer position for DCD frame counter progress meter.
- ![INFO][badge-info] Remove deprecated `which_types` function.

Version 2.14.1
-------------
- ![ENHANCEMENT][badge-enhancement] When `lastframe` is set and using `DCD` trajectory file, the initial trajectory reading will stop at `lastframe`. 
- ![BUGFIX][badge-bugfix] Fix bug in the final update of `coordination_number` when the random site count was zero (only appearing in very small trajectory tests).

Version 2.14.0
-------------
- ![FEATURE][badge-feature] Support for general 1-dimensional arrays as `frame_weights`. 

Version 2.13.2
-------------
- ![ENHANCEMENT][badge-enhancement] Improve progress meter of frame count in DCD files.

Version 2.13.1
-------------
- ![INFO][badge-info] Update python script, following the new selection features of PDBTools.jl 3.1.0.