# [Coordination numbers](@id Tools)

- [Computing radial distribution functions](@ref radial_distribution)
- [Overview of the solvent and solute properties](@ref overview)

```@meta
CollapsedDocStrings = true
```

## [Computing radial distribution functions](@id radial_distribution)

The distributions returned by the `mddf` function (the `mddf` and
`rdf` vectors), are normalized by the random reference state or
using a site count based on the numerical integration of the volume
corresponding to each minimum-distance to the solute. 

If, however, the
solute is defined by a single atom (as the oxygen atom of water, for
example), the numerical integration of the volume can be replaced by a
simple analytical spherical shell volume, reducing noise. The `ComplexMixtures.gr`
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
R = mddf(trajectory_file, solute, solvent, options)
g, kb = ComplexMixtures.gr(R.d,R.rdf_count,R.density.solvent_bulk,R.options.binstep)
```

```@autodocs
Modules = [ComplexMixtures]
Pages = ["gr.jl"]
```

## [Overview of the solvent and solute properties](@id overview)

The output to the REPL of the Result structure provides an overview of the
properties of the solution. The data can be retrieved into a data structure
using the `overview` function. Examples:     

```julia-repl
...

julia> results = mddf(trajectory_file, solute, solvent, Options(bulk_range=(8.0, 12.0)))

julia> results
--------------------------------------------------------------------------------
MDDF Overview - ComplexMixtures - Version 2.0.8
--------------------------------------------------------------------------------

Solvent properties:
-------------------

Simulation concentration: 0.49837225882780106 mol L⁻¹
Molar volume: 2006.532230249041 cm³ mol⁻¹

Concentration in bulk: 0.5182380507741433 mol L⁻¹
Molar volume in bulk: 1929.6151614228274 cm³ mol⁻¹

Solute properties:
------------------

Simulation Concentration: 0.002753437894076249 mol L⁻¹
Estimated solute partial molar volume: 13921.98945754469 cm³ mol⁻¹

Bulk range: 8.0 - 12.0 Å
Molar volume of the solute domain: 34753.1382279134 cm³ mol⁻¹

Auto-correlation: false

Trajectory files and weights:

   /home/user/NAMD/trajectory.dcd - w = 1.0

Long range MDDF mean (expected 1.0): 1.0378896753018338 ± 1.0920172247127446
Long range RDF mean (expected 1.0): 1.2147429551790854 ± 1.2081838161780682

--------------------------------------------------------------------------------
```

In this case, since solute and solvent are equivalent and the system is
homogeneous, the molar volumes and concentrations are similar. This is
not the case if the molecules are different or if the solute is at
infinite dilution (in which case the bulk solvent density might be
different from the solvent density in the simulation). 

To retrieve the data of the overview structure use, for example:

```julia-repl
julia> overview = overview(results);

julia> overview.solute_molar_volume
657.5051512801567
```
