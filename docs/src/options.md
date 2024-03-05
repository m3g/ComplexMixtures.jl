# [Options](@id options)

There are some options to control what exactly is going to be computed
to obtain the MDDF. These options can be defined by the user and passed to the
`mddf` function, using, for example: 

```julia
options = Options(lastframe=1000, bulk_range=(8.0, 12.0))
results = mddf(trajectory, options)
```

## Frame ranges and histogram properties

These are common options that the regular user might want to set 
in their calculation.

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

## Lower level options

These will probably never be set by the user, unless if dealing with 
some special system (large very large, or very low density system).

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

`seed`: Seed for random number generator. If `-1`, the seed will be
generated from the entropy of the system. If your results are dependent
on the seed, is is probable that you do not have enough sampling. Mostly
used for testing purposes. Two runs are only identical if ran with
the same seed and in serial mode.   

`StableRNG` (`::Bool`), defaults to `false`. Use a stable random number
generator from the `StableRNGs` package, to produce identical runs on
different architectures and Julia versions. Only used for testing. 

`nthreads`: How many threads to use. By default, it will be the number
of physical cores of the computer.
 
`lcell`: Integer, the cell length of the linked-cell method (actually
the cell length is `cutoff/lcell`). Default: 1.  

`GC`: Bool, force garbage collection, to avoid memory
overflow. Default: `true`. That this might be required is probably a result of
something that can vastly improved in memory management. This may slow down
parallel runs significantly if the GC runs too often.

`GC_threshold`: Float64, minimum fraction of the total memory of the
system required to force a GC run. That is, if `GC_threshold=0.1`, which
is the default, every time the free memory becomes less or equal to 10%
of the total memory available, a GC run occurs.  

## Frame statistical reweighing  

!!! compat
    Frame reweighing is available in ComplexMixtures 2.0.0 or greater.

Most times the weights of each frame of the trajectory are the same, resulting
from some standard MD simulation. If, for some reason, the frames have 
different statistical weights, the weights can be passed to the as an 
optional parameter `frame_weights`.

For example:
```julia-repl
julia> results = mddf(trajectory, options; frame_weights=[0.0, 1.0, 2.0])
```
The code above will assign a larger weight to the third frame of the trajectory.
These weights are relative (meaning that `[0.0, 1.0, 2.0]` would produce 
the same result). What will happen under the hood is that the distance counts
of the frames will be multiplied by each frame weight, and normalized for the
sum of the weights.

**Important:** The length of the `frame_weights` vector must be at least equal
to the number of the last frame read from the trajectory. That is, if `lastframe` 
is not set, and all the frames will be read, the length of `frame_weights` must
be equal to the length of the trajectory (the `stride` parameter will skip the
information both of the frames and its weights). If `lastframe` is set, then
the length of `frame_weights` must be at least `lastframe` (it can be greater,
and further values will be ignored). Importantly, the indices of the elements
in `frame_weights` are assumed to correspond to the indices of the frames
in the original trajectory file.

## Compute coordination number only

For some systems, it may be impossible, or to expensive, to compute the normalization
of the minimum-distance distribution function. Nevertheless, the coordination
number may still be an interesting information to be retrieved from the 
simulations. To run the computation to compute coordination numbers only, do:

```julia-repl
julia> results = mddf(trajectory, options; coordination_number_only = true)
```

!!! note    
    With `coordination_number_only` set to `true`, the arrays associated to
    MDDFs and KB integrals will be empty in the output data structure. 

```@autodocs
Modules = [ComplexMixtures]
Pages = ["Options.jl"]
```



