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

Version 2.19.1-DEV
--------------

Version 2.19.0
--------------
- ![FEATURE][badge-feature] `kbi(R; correction)`: running KBIs computed with the improved estimators of the infinite-volume KBI of Krüger and Vlugt (`correction=:G1` or `:G2`), which converge much faster than the truncated integral. Valid for radial distribution functions (single-atom solutes).
- ![FEATURE][badge-feature] `finite_volume_kbi` and `extrapolate_kbi`: finite-volume KBIs of spheres of diameter `L`, and their extrapolation to the infinite-volume limit, by fitting `G(L) = G∞ + F∞/L`.
- ![INFO][badge-info] New documentation page about the convergence and finite-size corrections of KBIs.
- ![FEATURE][badge-feature] Use PDBtools 3.41.0 visualization functions to document and display contributions interactivelly over structures.
- ![FEATURE][badge-feature] `volumetric_data`: converts the grid of `grid3D` into volumetric data (`PDBTools.VolumetricData`), optionally smoothed, which can be displayed as isosurfaces with `PDBTools.visualize` or written to OpenDX (`.dx`) files.
- ![FEATURE][badge-feature] `grid3D` supports solutes with multiple molecules: the grid is built around one of the molecules (keyword `molecule`, the first one by default).
- ![INFO][badge-info] Skip `:kbi` contribution computation if not available (because the read json is from an old version).
- ![INFO][badge-info] Use JSON.jl (v1) instead of the deprecated JSON3.jl (and StructTypes.jl) to read and write results files. The file format is unchanged.
- ![FIX][badge-fix] Fix the display of `TrajectoryFileOptions` (e.g. `results.files[1]`) with non-empty frame weights.
- ![INFO][badge-info] Update references and application papers.
- ![BUGFIX][badge-bugfix] `gr(R::Result)` now returns `R.rdf` and `R.kb_rdf` (normalized by the random reference state) when the solute has more than one atom per molecule. Previously, the spherical shell volume was used, which is only valid for single-atom solutes, and produced wrong g(r) and KB integrals in that case. The behavior for single-atom solutes is unchanged. 
- ![BUGFIX][badge-bugfix] Random rotations of solvent molecules in the ideal-gas reference state are now uniformly distributed. The previous sampling of the rotation angles was not uniform in orientation space. The practical relevance is small: the molecules are copied from the bulk before being rotated, so in isotropic solutions their orientations were already random and the bias cancelled. The bias could affect, slightly, the MDDF and KBI of multi-atom solvents (not the site-based `rdf`) if the bulk orientations are anisotropic or poorly sampled (e.g. few solvent molecules, short trajectories). Results for a given random seed will differ slightly from previous versions.

Version 2.18.2
--------------
- ![INFO][badge-info] Import `save` and `load` from MolSimToolkitShared (v1.6) for interop with PDBTools.jl. Note: this required the "benign" piracy of `save` from `MolSimToolkitShared`, to avoid breaking `save(::String)`. This will be fixed in a future v3.

Version 2.18.1
--------------
- ![INFO][badge-info] Fix concurrency issue of RNG generator and unnccessary lock of system update.

Version 2.18.0
--------------
- ![FEATURE][badge-feature] Provide the `type=:kbi` option to `contributions` and `ResidueContributions` functions, to obtain the proximal contributions to the KBIs. 
- ![FEATURE][badge-feature] Save `solute_group_count_random` and `solvent_group_count_random` to allow the decomposition of KBIs into proximal group contributions. 

Version 2.17.1
--------------
- ![INFO][badge-info] Update internals to use CellListMap v0.10.1. 

Version 2.17.0
--------------
- ![FEATURE][badge-feature] Support for `AbstractString` in the signature of functions that previously only accepted `String`.
- ![INFO][badge-info] add concepts.md page - use top_menu
- ![INFO][badge-info] Add JACS Au cover and update reference.

Version 2.16.0
--------------
- ![FEATURE][badge-feature] The `grid3D` function now supports the `type` keyword argument, to choose from `:mddf` (default), `:coordination_number`, or `:md_count`, parameter that is passed to `contributions` to build the grid.

Version 2.15.1
--------------
- ![ENHANCEMENT][badge-enhancement] Provide some better error messages for group contribution retrieval argument errors.
- ![INFO][badge-info] Update version of ShowMethodTesting and adjust tests.
- ![INFO][badge-info] Update PDBTools compat requirement to 3.
- ![INFO][badge-info] Remove testing temporary files.
- ![INFO][badge-info] Store test files in Zenodo server instead of Dropbox.

Version 2.15.0
--------------
- ![FEATURE][badge-feature] Implement the `renormalize` function to compute the normalization with different bulk densities, as using the bulk density as the total simulation density, or a custom density.
- ![INFO][badge-info] Update publication list.

Version 2.14.4
--------------
- ![ENHANCEMENT][badge-enhancement] Read chemfiles trajectories in place.

Version 2.14.3
--------------
- ![BUGFIX][badge-bugfix] fix diagonal unitcell test for when cell is not orthorhombic and contains negative vector entries.
- ![INFO][badge-info] Update FortranFiles dependency to 0.6.2 and use `seekstart` instead of `rewind`. 

Version 2.14.2
--------------
- ![ENHANCEMENT][badge-enhancement] Better finalizer position for DCD frame counter progress meter.
- ![INFO][badge-info] Remove deprecated `which_types` function.

Version 2.14.1
--------------
- ![ENHANCEMENT][badge-enhancement] When `lastframe` is set and using `DCD` trajectory file, the initial trajectory reading will stop at `lastframe`. 
- ![BUGFIX][badge-bugfix] Fix bug in the final update of `coordination_number` when the random site count was zero (only appearing in very small trajectory tests).

Version 2.14.0
--------------
- ![FEATURE][badge-feature] Support for general 1-dimensional arrays as `frame_weights`. 

Version 2.13.2
--------------
- ![ENHANCEMENT][badge-enhancement] Improve progress meter of frame count in DCD files.

Version 2.13.1
--------------
- ![INFO][badge-info] Update python script, following the new selection features of PDBTools.jl 3.1.0.