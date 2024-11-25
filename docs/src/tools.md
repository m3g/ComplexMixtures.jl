# [Tools](@id Tools)

- [Coordination numbers](@ref coordination_number)
- [2D density map per residue](@ref 2D_per_residue)
- [3D density map around a macromolecule](@ref grid3D)
- [Overview of the solvent and solute properties](@ref overview)
- [Computing radial distribution functions](@ref radial_distribution)

```@meta
CollapsedDocStrings = true
```

## [Coordination numbers](@id coordination_number)

The function
```julia
coordination_number(R::Result, group::Union{SoluteGroup, SolventGroup})
```
computes the coordination number of a given group of atoms from the solute or solvent atomic contributions to the MDDF. Here, `R` is
the result of the `mddf` calculation, and `group_contributions` is the output of the `contributions` function for the desired set of atoms.

If no group is defined, the coordination number of the complete solute is returned, which is equivalent to the `R.coordination_number` field of the `Result` data structure:
```
coordination_number(R::Result) == R.coordination_number
```

!!! note
    There are some systems for which the normalization of the distributions is not 
    necessary or possible. It is still possible to compute the coordination numbers,
    by running, instead of `mddf`, the `coordination_number` function:
    ```
    coordination_number(trajectory::Trajectory, options::Options)
    ```
    This call will return `Result` data structure but with all fields requiring 
    normalization with zeros. In summary, this result
    data structure can be used to compute the coordination numbers, but not the MDDF, RDF, or KB integrals.

!!! compat
    The independent computation of coordination numbers was introduced in version 1.1.

### Example

In the following example we compute the coordination number of the atoms of `residue 50` (which belongs to the solute - a protein) with the solvent atoms of TMAO, as a function of the distance. The plot produced will show side by side the residue contribution to the MDDF and the corresponding coordination number.

```julia
using ComplexMixtures, PDBTools
using Plots, EasyFit
pdb = readPDB("test/data/NAMD/structure.pdb")
R = load("test/data/NAMD/protein_tmao.json")
solute = AtomSelection(PDBTools.select(pdb, "protein"), nmols=1)
residue50 = PDBTools.select(pdb, "residue 50")
# Compute the group contribution to the MDDF
residue50_contribution = contributions(R, SoluteGroup(residue50))
# Now compute the coordination number
residue50_coordination = coordination_number(R, SoluteGroup(residue50))
# Plot with twin y-axis
plot(R.d, movavg(residue50_contribution,n=10).x,
    xaxis="distance / Å", 
    yaxis="MDDF contribution", 
    linewidth=2, label=nothing, color=1
)
plot!(twinx(),R.d, residue50_coordination, 
    yaxis="Coordination number", 
    linewidth=2, label=nothing, color=2
)
plot!(title="Residue 50", framestyle=:box, subplot=1)
```

With appropriate input data, this code produces:

```@raw html
<center>
<img width=60% src="../figures/coordination.png" width=80%>
</center>
```

```@autodocs
Modules = [ComplexMixtures]
Pages = ["coordination_number.jl"]
```

## [2D density map per residue](@id 2D_per_residue)

- [The `ResidueContributions` object](@ref)
- [Indexing, slicing, arithmetic operations](@ref)
- [Saving and loading a ResidueContributions object](@ref)

### The `ResidueContributions` object

One nice way to visualize the accumulation or depletion of a solvent around a macromolecule (a protein, for example), is to obtain a 2D map of the density as a function of the distance from its surface. For example, in the figure below the density of a solute (here, Glycerol), in the neighborhood of a protein is shown:

```@raw html
<center>
<img src="../figures/density.png" width=80%>
</center>
```

Here, one can see that Glycerol accumulates on Asp76 and on the proximity of hydrogen-bonding residues (Serine residues mostly). This figure was obtained by extracting from atomic contributions of the protein the contribution of each residue to the MDDF, coordination numbers or minimum-distance counts. 

!!! compat
    All features described in this section are only available in v2.8.0 or greater.

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

The `ResidueContribution` object can be used to produce a high-quality contour plot using the `Plots.contourf` function:

```@docs
Plots.contourf(::ResidueContributions)
```
A complete example of its usage can be seen [here](@ref 2D-map-example1). 

### Indexing, slicing, arithmetic operations

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

### Saving and loading a ResidueContributions object

The `ResidueContributions` object can be saved and loaded for easier data analysis. In particular, this 
is important for very large structures, where its computation can be costly. The saving and loading 
functions can be use with:

```julia
rc = ResidueContributions(results1, select(atoms, "protein")); 
# Save rc objecto to a file (json format):
save("residue_contributions.json", rc) 
# Load json file into a new rc_loaded object:
rc_loaded = load("residue_contributions.json", ResidueContributions)
```

Note that the `load` function requires, as a second argument, the `ResidueContributions` type, to differentiate
the method from the loading of the `Result` data structure.

```@docs
load(::String, ::Type{ResidueContributions})
save(::String, ::ResidueContributions)
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
    
## [3D density map around a macromolecule](@id grid3D)

Three-dimensional representations of the distribution functions can also be obtained from the MDDF results. These 3D representations are obtained from the fact that the MDDFs can be decomposed into the contributions of each solute atom, and that each point in space is closest to a single solute atom as well. Thus, each point in space can be associated to one solute atom, and the contribution of that atom to the MDDF at the corresponding distance can be obtained.   

A 3D density map is constructed with the `grid3D` function:

```@autodocs
Modules = [ComplexMixtures]
Pages = ["tools/grid3D.jl"]
```

The call to `grid3D` will write an output a PDB file with the grid points, which loaded in a visualization software side-by-side with the protein structure, allows the production of the images shown. The `grid.pdb` file contains a regular PDB format where: 

- The positions of the atoms are grid points. 
- The identity of the atoms correspond to the identity of the protein atom contributing to the MDDF at that point (the closest protein atom). 
- The temperature-factor column (`beta`) contains the relative contribution of that atom to the MDDF at the corresponding distance. 
- The `occupancy` field contains the distance itself.

For example, the distribution function of a hydrogen-bonding liquid solvating a protein will display a characteristic peak at about 1.8Å. The MDDF at that distance can be decomposed into the contributions of all atoms of the protein which were found to form hydrogen bonds to the solvent. A 3D representation of these contributions can be obtained by computing, around a static protein (solute) structure, which are the regions in space which are closer to each atom of the protein. The position in space is then marked with the atom of the protein to which that region "belongs" and with the contribution of that atom to the MDDF at each distance within that region. A special function to compute this 3D distribution is provided here: `grid3D`. 

This is better illustrated by a graphical representation. In the figure below we see a 3D representation of the MDDF of Glycerol around a protein, computed from a simulation of this protein in a mixture of water and Glycerol. A complete set of files and a script to reproduce this example [is available here](@ref 3D-map-example1). 

```@raw html
<center>
<img src="../figures/density3D_final.png" width=100%>
</center>
```

In the figure on the left, the points in space around the protein are selected with the following properties: distance from the protein smaller than 2.0Å and relative contribution to the MDDF at the corresponding distance of at least 10% of the maximum contribution. Thus, we are selecting the regions of the protein corresponding to the most stable hydrogen-bonding interactions. The color of the points is the contribution to the MDDF, from blue to red. Thus, the most reddish-points corresponds to the regions where the most stable hydrogen bonds were formed. We have marked two regions here, on opposite sides of the protein, with arrows.

Clicking on those points we obtain which are the atoms of the protein contributing to the MDDF at that region. In particular, the arrow on the right points to the strongest red region, which corresponds to an Aspartic acid. These residues are shown explicitly under the density (represented as a transparent surface) on the figure in the center.   

The figure on the right displays, overlapped with the hydrogen-bonding residues, the most important contributions to the second peak of the distribution, corresponding to distances from the protein between 2.0 and 3.5Å. Notably, the regions involved are different from the ones forming hydrogen bonds, indicating that non-specific interactions with the protein (and not a second solvation shell) are responsible for the second peak. 

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

## Legacy 

The legacy function `contourf_per_residue` provides a direct way to produce contour plots:

```@docs
contourf_per_residue(::Result, ::AbstractVector{PDBTools.Atom})
```


