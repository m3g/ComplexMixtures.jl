```@meta
CollapsedDocStrings = true
```

# [2D density maps](@id 2D_per_residue)

- [The `ResidueContributions` object](@ref)
- [Contributions of subgroups of residues](@ref)
- [Indexing, slicing, arithmetic operations](@ref)
- [Saving and loading a ResidueContributions object](@ref)

## The `ResidueContributions` object

One nice way to visualize the accumulation or depletion of a solvent around a macromolecule (a protein, for example), is to obtain a 2D map of the density as a function of the distance from its surface. For example, in the figure below the density of a solute (here, Glycerol), in the neighborhood of a protein is shown:

```@raw html
<center>
<img src="../figures/density.png" width=80%>
</center>
```

Here, one can see that Glycerol accumulates on Asp76 and on the proximity of hydrogen-bonding residues (Serine residues mostly). This figure was obtained by extracting from atomic contributions of the protein the contribution of each residue to the MDDF, coordination numbers or minimum-distance counts. 

The computation of the contributions of each residue can be performed with the convenience function `ResidueContributions`, which
creates an object containing the contributions of the residues to the mddf (or coordination numbers, or minimum-distance counts), the 
residue names, and distances:

```@docs
ResidueContributions
```

The output of `ResidueContributions` is by default shown as a simple unicode plot:

```@raw html
<center>
<img src="../figures/ResidueContributions.png" width=60%>
</center>
```

The `ResidueContribution` object can be used to produce a high-quality contour plot using the `Plots.contourf` (or `contour`, or `heatmap`) function:

```@docs
Plots.contourf(::ResidueContributions)
```

A complete example of its usage can be seen [here](@ref 2D-map-example1). 

!!! compat
    All features described in this section are only available in v2.10.0 or greater.

## Contributions of subgroups of residues

Residue contributions can also be computed for subgroups of the residues. For example, as a continuation of the [protein in glycerol example](@ref 2D-map-example1) , one
can compute the map of residue contributions, but splitting the contributions of backbone and side-chains of the residues:

```julia
rc_bb = ResidueContributions(
  results, 
  select(protein, "backbone and resnum >= 70 and resnum <= 110")
)
rc_st = ResidueContributions(
  results, 
  select(protein, "sidechain and resnum >= 70 and resnum <= 110")
)
```

And we plot the difference between these two maps:
```julia
using Plots
contourf(rc_st - rc_bb; oneletter=true)
```

obtaining the following figure:

```@raw html
<center>
<img width=70% src="../assets/scripts/example1/2D-map_st-bb.png">
</center>
```
which shows that the side-chains contribute mostly to these densities, except, as expected, for some Gly residues.

## Indexing, slicing, arithmetic operations

The `ResidueContributions` object can be indexes and sliced, for the analysis of the contributions of specific residues
or range of residues:

```julia
rc = ResidueContributions(results1, select(atoms, "protein")); 
rc_7 = rc[7] # contributions of residue 7
rc_range = rc[20:50] # contributions of a range of residues
```
Slicing will return a new `ResidueContributions` object.

Additionally, these `ResidueContributions` objects can be subtracted, divided, summed, or multiplied, to compare contributions of residues
among different simulations. Typically, if one wants to compare the solvation of residues in two different simulations, 
one can do:
```julia
# first simulation (for example, low temperature)
rc1 = ResidueContributions(results1, select(atoms, "protein")); 

# second simulation (for example, high temperature)
rc2 = ResidueContributions(results2, select(atoms, "protein"));

# difference in residue contributions to solvation
rc_diff = rc2 - rc1

# Plot difference
using Plots
contourf(rc_diff; title="Density difference", step=2, colorbar=:left)
```
Which will produce a plot similar to the one below (the data of this plot is just illustrative):

which will return a new `ResidueContributions` object.

```@raw html
<center>
<img src="../figures/density2.png" width=70%>
</center>
```
Finally, it is also possible to renormalize the contributions by multiplication or division by scalars,
```julia
rc2 = rc / 15
rc2 = 2 * rc
```

When the contributions of a single residue are computed, or a single-residue contribution is retrieved from
a `ResidueContributions` object, the indexing and iteration over that object occurs over the contributions of that residue:

```julia
using ComplexMixtures, PDBTools, Plots
...
result = mddf(trajectory_file, solute, solvent, options)
rc = ResidueContributions(result, select(atoms, "protein"))
rc7 = rc[7] # contributions of residue 7
# iterate over the contributions of residue 7
rc7[1] # contribution of the first distance
rc7[end] # contribution of the last distance
```

This is particular useful to retrieve the contributions from all residues at a given distance:

```julia
rc = ResidueContributions(result, select(atoms, "protein"))
rc_last_distance = [ r[end] for r in rc ] 
# or, equivalently
rc_last_distance = last.(rc)
# compute the maximum contribution of each residue:
max_c = maximum.(rc)
```

## Saving and loading a ResidueContributions object

The `ResidueContributions` object can be saved and loaded for easier data analysis. In particular, this 
is important for very large structures, where its computation can be costly. The saving and loading 
functions can be use with:

```julia
rc = ResidueContributions(results1, select(atoms, "protein")); 
# Save rc object to a file (json format):
save("residue_contributions.json", rc) 
# Load json file into a new rc_loaded object:
rc_loaded = load("residue_contributions.json", ResidueContributions)
```

Note that the `load` function requires, as a second argument, the `ResidueContributions` type, to differentiate
the method from the loading of the `Result` data structure.

```@docs
load(::AbstractString, ::Type{ResidueContributions})
save(::AbstractString, ::ResidueContributions)
```

!!! tip 
    These `ResidueContributions` methods are convenience functions only. 

    Basically, we are extracting the contribution of each residue independently and building a matrix where each row 
    represents a distance and each column a residue.  Using `PDBTools`, this can be done with, for example: 
    
    ```julia
    residues = collect(eachresidue(protein))
    residue_contributions = zeros(length(R.d),length(residues))
    for (i,residue) in pairs(residues)
      c = contributions(results, SoluteGroup(residue)) 
      residue_contributions[:,i] .= c
    end
    ```
    
    The above produces a matrix with a number of columns equal to the number of residues and a number of rows equal to the number of MDDF points. That matrix can be plotted as a contour map with adequate plotting software. 
    
