# [Tools](@id Tools)

Some tools are provided to analyze the results:

## Overview of the solvent and solute properties 

The output to the REPL of the Result structure provides an overview of the
properties of the solution. The data can be retrieved into a data structure
using the `overview` function. Examples:     

```julia
...
julia> results = ComplexMixtures.mddf(trajectory)

julia> results

-------------------------------------------------------------------------------

 MDDF Overview: 

 Solvent properties: 
 ------------------- 

 Simulation concentration: 1.5209006318095133 mol L⁻¹
 Molar volume: 657.5051512801567 cm³ mol⁻¹

 Concentration in bulk: 1.4918842545752287 mol L⁻¹
 Molar volume in bulk: 670.2932864484995 cm³ mol⁻¹ 

 Solute properties: 
 ------------------ 

 Simulation Concentration: 1.5209006318095133 mol L⁻¹
 Estimated solute partial molar volume: 657.5051512801567 cm³ mol⁻¹

 Using with dbulk = 20.0Å: 
 Molar volume of the solute domain: 30292.570006549242 cm³ mol⁻¹

 Auto-correlation: true

 Trajectory files and weights: 
   ./vinicius.xtc - w = 1.0

 Long range MDDF mean (expected 1.0): 1.1090804621839963 +/- 0.04298849642932878
 Long range RDF mean (expected 1.0): 1.15912932236198 +/- 0.05735018864444404

-------------------------------------------------------------------------------

```

In this case, since solute and solvent are equivalent and the system is
homogeneous, the molar volumes and concentrations are similar. This is
not the case if the molecules are different or if the solute is at
infinite dilution (in which case the bulk solvent density might be
different from the solvent density in the simulation). 

To retrieve the data of the overview strcture use, for example:

```julia
julia> overview = ComplexMixtures.overview(results);

julia> overview.solute_molar_volume
657.5051512801567

```

## Computing a 2D density map around a macromolecule 

One nice way to visualize the accumulation or depletion of a solvent around a macromolecule (a protein, for example), is to obtain a 2D map of the density as a function of the distance from its surface. For example, in the figure below the density of a solute (here, Glycerol), in the neighborhood of a protein is shown:

```@raw html
<center>
<img src="../figures/density.png" width=80%>
</center>
```

Here, one can see that Glycerol accumulates on Asp76 and on the proximity of hydrogen-bonding residues (Serine residues mostly). This figure was obtained by extracting from atomic contributions of the protein the contribution of each residue to the MDDF. Using `PDBTools`, this can be done with, for example: 

```julia
residues = collect(eachresidue(protein))
residue_contributions = zeros(length(R.d),length(residues))
for (i,residue) in pairs(residues)
  c = ComplexMixtures.contrib(solute,R.solute_atom,residue) 
  residue_contributions[i,:] .= c
end
```

The above produces a matrix with a number of columns equal to the number of residues and a number of rows equal to the number of MDDF points. That matrix can be plotted as a contour map with adequate plotting software. [A complete running example is provided here](https://github.com/m3g/ComplexMixturesExamples/tree/main/Density2D), producing the figure above.    

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




