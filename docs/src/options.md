# [Options](@id options)

There are some options to control what exactly is going to be computed
to obtain the MDDF. These options can be defined by the user and passed to the
`mddf` function, using, for example: 

```julia
options = ComplexMixtures.Options(lastframe=1000)
results = ComplexMixtures.mddf(trajectory,options)
```

### Most common options that the user might want to control are:

`firstframe`: Integer, first frame of the trajectory to be considered.

`lastframe`: Integer, last frame of the trajectory to be considered.

`stride`: Integer, consider every stride frames, that is, if `stride=5`
only one in five frames will be considered.

`binstep`: Real, length of the bin step of the histograms, default =
0.02 Angstroms.

`dbulk`: Real, distance from which the solution is to be considered as a
bulk solution, that is, where the presence of the solute does not affect
the structure of the solution anymore. This parameter is important in
particular for systems with a single solute molecule (a protein, for
example), where the density of the solvent in the box is not the bulk
density of the solvent, which must be computed independently. Default:
10 Angstroms. 

`cutoff`: Real, the maximum distance to be considered in the
construction of histograms. Default: 10 Angstroms. 

`usecutoff`: `true/false`: If true, the cutoff distance might be
different from `dbulk` and the density of the solvent in bulk will be
estimated from the density within `dbulk` and `cutoff`. If `false`, the
density of the solvent is estimated from the density outside `dbulk` by
exclusion. Default: `false`. 

### Options that most users will probably never change:

`irefatom`: Integer, index of the reference atom in the solvent molecule
used to compute the shell volumes and domain volumes in the Monte-Carlo
volume estimates. The final `rdf` data is reported for this atom as
well. By default, we choose the atom which is closer to the center of
coordinates of the molecule, but any choice should be fine. 

`n_random_samples`: Integer, how many samples of random molecules are
generated for each solvent molecule to compute the shell volumes and
random MDDF counts. Default: 10. Increase this only if you have short
trajectory and want to obtain reproducible results for that short
trajectory. For long trajectories (most desirable and common), this
value can even be decreased to speed up the calculations. 

`lcell`: Integer, the cell length of the linked-cell method (actually
the cell length is `cutoff/lcell`). Default: 2.  

`sleep`: Real, the time between checks between multiple spawns of
calculations in parallel. Default 0.1 s. 

`GC`: Bool, force garbage collection in parallel runs, to avoid memory
overflow. Default: `true`. That this might be required is probably a result of
something that can vastly improved in memory management. This may slow down
parallel runs significantly if the GC runs too often.

`GC_threshold`: Float64, minimum fraction of the total memory of the
system required to force a GC run. That is, if `GC_threshold=0.1`, which
is the default, every time the free memory becomes less or equal to 10%
of the total memory available, a GC run occurs.  

