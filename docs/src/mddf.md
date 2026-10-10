```@meta
CollapsedDocStrings = true
```
# Computing the MDDF

## Minimum-Distance Distribution Function

The main function of the ComplexMixtures package actually computes the MDDF between
the solute and the solvent chosen. 

```@docs
mddf
```

The `mddf` functions is run with, for example:

```julia
results = mddf(trajectory_file, solute, solvent, Options(cutoff=15.0))  
```

The MDDF along with other results, like the corresponding KB integrals,
are returned in the `results` data structure, which is described in the
[next section](@ref results).

It is possible to tune several options of the calculation, by setting
the `Options` data structure with user-defined values in advance.
The most common parameters to be set by the user are `cutoff`
and `stride`. 

`stride` defines if some frames will be skip during the calculation (for
speedup). For example, if `stride=5`, only one in five frames will be
considered. Adjust stride with:  

```julia
options = Options(stride=5, cutoff=15.0)
results = mddf(trajectory_file, solute, solvent, options)
```

!!! note
    `cutoff` defines the maximum distance to the solute for which the distribution
    function is computed. The bulk density of the solvent is estimated from the region of
    the system beyond the cutoff. Thus, the cutoff must be large enough such that,
    beyond it, the solute does not significantly affect the structure of the solvent anymore. 
    The adequate choice of `cutoff` can be inspected by the convergence of the distribution 
    functions (which must converge to 1.0), by the convergence of the KB integrals, 
    and with the [`reference_density`](@ref) function. By default, `cutoff=10.0`, but it is 
    *highly recommended* to set it according to the system.

    The `bulk_range` option, which defined a range of distances from which the bulk density was 
    estimated, is deprecated since version 2.19.1.

See the [Options](@ref options) section for further details and other options
to set.

## Coordination numbers only

The coordination number is the unnormalized count of how many molecules of the solvent
are within a given distance to the solute. Coordination numbers can be computed 
for systems where the normalization of the distribution functions is not possible
(or needed) because of an ill definition of an ideal-gas state. For example, 
in highly heterogeneous systems, on in systems with only a few molecules of the 
"solvent", the density of the bulk solution might not be properly defined.   

In these cases, nevertheless, coordination numbers can be computed and still 
provide valuable information about the molecular structure of the system. Coordination
number can be computed also from the results obtained from a `mddf` run, as explained in 
the [corresponding section of the Tools menu](@ref coordination_number).

The `coordination_number` function, called with the same arguments as the `mddf`
function, can be used to compute coordination numbers without the normalization
required for the MDDF, providing (possibly much) faster computations when the 
normalization is not possible or required:

```@docs
coordination_number(::AbstractString, ::AtomSelection, ::AtomSelection, ::Options)
```

!!! note 
    The `mddf`, `kb`, and random count parameters will be filled with zeros when using 
    this options, and are meaningless. 
