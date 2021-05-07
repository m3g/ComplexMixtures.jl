# Computing the Minimum-Distance Distribution Function

The main function of the ComplexMixtures package actually computes the MDDF between
the solute and the solvent chosen. 

It is run with the following command:

```julia
results = mddf(trajectory)  
```

The MDDF along with other results, like the corresponding KB integrals,
are returned in the `results` data structure, which is described in the
[next section](@ref results).

It is possible to tune several options of the calculation, by setting
the `Options` data structure with user-defined values in advance.
The most common parameters to be set by the user are probably `dbulk`
and `stride`. 

`dbulk` defines the distance from the solute above which
the user believes that the reference solute molecule does not
significantly anymore the structure of the solvent. The default value is
10 Angstroms, but for large solvent molecules this might not be enough.
To increase dbulk, use:  
```julia
options = Options(dbulk=15.)
results = mddf(trajectory,options)
```

`stride` defines if some frames will be skip during the calculation (for
speedup). For example, if `stride=5`, only one in five frames will be
considered. Adjust stride with:  
```julia
options = Options(stride=5)
results = mddf(trajectory,options)
```

See the [Options](@ref options) section for further details and other options
to set.

