# Computing the Minimum-Distance Distribution Function

The main function of the ComplexMixtures package actually computes the MDDF between
the solute and the solvent chosen. 

It is run with the following command:

```julia
results = mddf(trajectory, Options(bulk_range=(10.0, 15.0)))  
```

The MDDF along with other results, like the corresponding KB integrals,
are returned in the `results` data structure, which is described in the
[next section](@ref results).

It is possible to tune several options of the calculation, by setting
the `Options` data structure with user-defined values in advance.
The most common parameters to be set by the user are `bulk_range`
and `stride`. 

`stride` defines if some frames will be skip during the calculation (for
speedup). For example, if `stride=5`, only one in five frames will be
considered. Adjust stride with:  

```julia
options = Options(stride=5, bulk_range=(10.0, 15.0))
results = mddf(trajectory, options)
```

!!! note
    `bulk_range` defines the subset of the system, as defined according
    to a range of distances from the solute, that are to be considered
    as the bulk solution. Within this range of distances, the user 
    believes that the reference solute molecule does not
    significantly affect anymore the structure of the solvent. 

    By default, all molecules above 10 Angstroms from the solute are
    considered bulk molecules (corresponding to `Options(dbulk=10.0)`), but
    it is *highly recommended* to use a manual definition of `bulk_range`.

    The definition of a range of distances within the system to compute the
    bulk density is adequate because this system subset is then an open
    system with a solvent molecule reservoir. The adequate choice of `bulk_range`
    can be inspected by the proper convergence of the distribution functions
    (which must converge to 1.0) and a proper convergence of the KB integrals.

    The `bulk_range` option was introduced in version 2.1.0.

See the [Options](@ref options) section for further details and other options
to set.

## Reference functions

```@autodocs
Modules = [ComplexMixtures]
Pages = ["mddf.jl"]
```
