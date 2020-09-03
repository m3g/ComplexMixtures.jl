# [Tools](@id Tools)

Some tools are provided to analyze the results:

## Overview of the solvent and solute properties 

The overview function displays global properties of the solute and
solvent molecules in the simulation, given the results. For example,
this is the output from a simple simulation of pure water:

```julia
julia> ComplexMixtures.overview(R)

 Overview: 

 Solvent properties: 
 ------------------- 

 Simulation concentration: 51.398872485722066 mol L⁻¹
 Molar volume: 19.455679699545684 cm³ mol⁻¹

 Concentration in bulk: 56.20970009706161 mol L⁻¹
 Molar volume in bulk: 17.790523668925882 cm³ mol⁻¹ 

 Solute properties: 
 ------------------ 

 Simulation Concentration: 56.13503056725589 mol L⁻¹
 Simulation molar volume: 17.814188215358516 cm³ mol⁻¹

 Using with dbulk=15.0 Å: 
 Volume of the solute domain: 8492.91003652388 cm³ mol⁻¹
 Molar volume of the solute domain: 5.114550052577582e54 mol L⁻¹

```

In this case, since solute and solvent are equivalent and the system is
homogeneous, the molar volumes and concentrations are similar. This is
not the case if the molecules are different or if the solute is at
infinite dilution (in which case the bulk solvent density might be
different from the solvent density in the simulation). 

## Computing radial distribution functions

The distributions returned by the `mddf` function (meaning `mddf` and
`rdf`) vectors, are normalized either the random reference state or
using a site count based on the numerical integration of the volume
corresponding to each minimum-distance to the solute. If, however, the
solute is defined by a single atom (as the oxygen atom of water, for
example), the numerical integration of the volume can be replaced by a
simple analytical spherical shell volume, reducing noise. The `gr`
function returns the radial distribution function and the KB integral 
computed from the results, using this volume estimate: 

```julia
g, kb = ComplexMixtures.gr(R)

```

By default, the single-reference count (`rdf_count`) of the Result
structure will be used to compute the radial distribution function. The
function can be called with explicit control of all input parameters: 

```julia
g, kb = ComplexMixtures.gr(r,count,density,binstep)

```
where:

| Parameter | Definition | Result structure output data to provide |
|:---------:|:-----------|:-------------------|
| `r`       | Vector of distances | The `d` vector |
| `count`   | Number of site counts at each r | The `rdf` or `mddf` vectors  |
| `density` | Bulk density | The `density.solvent_bulk` or `density.solvent` densities. |
| `binstep` | The histogram step | The `options.binstep`  |
|           |                    |                        |

Example:
```julia
...
R = ComplexMixtures.mddf(trajectory,options)
g, kb = ComplexMixtures.gr(R.d,R.rdf_count,R.density.solvent_bulk,R.options.binstep)

```




