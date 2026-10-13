# [Version compatibility](@id version-compatibility)

## [Migrating from 2.19 to 3.0](@id migration-3)

Scripts written for version 2.19 run with version 3.0, but some defaults changed, and the
numerical results change. The most important changes are:

### The number of random samples (`n_random_samples`)

!!! warning
    Do not use large values of `n_random_samples` (as was common in previous versions) for 
    self-correlations, or for cross-correlations with more than one solute molecule. In these 
    cases each random sample now costs about as much as the computation of the minimum distances 
    of the frame itself, such that, for example, `n_random_samples=300` makes the calculation 
    about two orders of magnitude slower than with the default (1). Remove the option from the
    `Options`, to use the default.

In version 3.0, each random (ideal-gas) sample of self-correlations (the solute and solvent are
the same molecules) and of cross-correlations with more than one solute molecule (for example,
glycerol in water/glycerol mixtures) is a configuration in which all molecules are placed at 
random, and the minimum distances to all the solute molecules are computed, as for the actual 
frame. Thus, one sample provides one reference for each solute molecule, and the default is 
`n_random_samples=1`. Previously, each sample provided a single reference molecule, such that 
the normalizations, and particularly the KB integrals, were noisy, which is why large numbers of
samples were often used.

For cross-correlations with a single solute molecule (for example, a protein), nothing changed: 
each sample provides one reference, and the default is 10 samples per frame.

See the [Options](@ref options) section for details. 

### Kirkwood-Buff integrals

- `kbi(R)` is the recommended way to obtain the KB integrals. Its defaults are now 
  `correction=:W7` (a weight function that damps the oscillations of the running integral at 
  long distances) and `normalization=:ganguly` (the reference density is the density of the 
  solvent outside the solute domain at each distance). `R.kb` is now equal to `kbi(R)`.
- The KB integrals of version 2.19 (`R.kb`, the truncated integral of the MDDF) are obtained with 
  `kbi(R; correction=:none, normalization=:mddf)`.
- `contributions(...; type=:kbi)` and `ResidueContributions(...; type=:kbi)` accept the same 
  `correction` and `normalization` options, and, by default, sum to `kbi(R)`.
- For results saved with previous versions, `R.kb` keeps the stored values: use `kbi(R)` to 
  compute the corrected KB integrals from them.

See the [Kirkwood-Buff integrals](@ref kbi) section for details.

### Bulk density and cutoff

- The bulk density of the solvent is always estimated from the region beyond the `cutoff`, 
  and the only option that defines it is the `cutoff`. The `bulk_range`, `dbulk`, and 
  `usecutoff` options are deprecated (a warning is issued), and mapped to the `cutoff`. 
  Results computed with `bulk_range` change, because the bulk density is estimated from a 
  different region.
- The solvent molecules are classified as bulk molecules by their minimum distance to the 
  solute (instead of by the distance of their reference atom). For solvents with more than 
  one atom, the bulk density, and thus `R.mddf`, change slightly.
- The new [`reference_density`](@ref) function can be used to inspect the uniformity of the
  density of the solvent beyond the cutoff (see [Inspecting the reference density](@ref kbi_reference_density)).

### Changes in the results

- The MDDFs of self-correlations were overestimated by a factor ``N/(N-1)``, where ``N`` is the
  number of molecules (for example, by 0.55% for 181 molecules). This was fixed.
- The results of self-correlations and of cross-correlations with more than one solute molecule
  change within their statistical uncertainty, because of the new random samples, which are 
  also more converged.

## [Migrating from 1.X to 2.X](@id updating-scripts)

The syntax changes necessary to update scripts from version `1.X` to `2.X` of 
the package are:

### Atom selections

The previous `Selection` structure was renamed to `AtomSelection` for clarity.
- Before:
```julia
water = Selection(water; natomspermol=3)
```
- Now:
```julia
water = AtomSelection(water; natomspermol=3)
```

### Group contributions syntax

The syntax to computing group contributions is improved. Previously, the `contrib` or
`contributions` functions required three somewhat redundant parameters. 
- Before:
The call to `contributions` required 3 parameters: the `Selection` structure,
the matrix of contributions, and the indexes of the atoms for which the
contributions were desired:
```julia
h_contributions = contributions(solvent, R.solvent_atom, h_indexes)
```
- Now:
The contributions are extracted from the `Result` data structure, by 
providing either a `SoluteGroup` or `SolventGroup` object, which are
setup with the group names, group indexes, atom names, or atom indexes:
```julia
h_contributions = contributions(R, SolventGroup(h_indexes))
```

### Frame weights

`frame_weights` is an option of the `mddf` execution. That is previously,
they were defined in the `Options` data structure, and now they are passed
to the `mddf` function.
- Before:
```julia
options = Options(frame_weights=[1.0, 2.0], bulk_range=(8.0, 12.0))
results = mddf(trajectory_file, solute, solvent, options)
```
- Now:
```julia
results = mddf(trajectory_file, solute, solvent, options; frame_weights=[1.0, 2.0])
```
