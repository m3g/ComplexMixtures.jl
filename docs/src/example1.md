# Protein in water/glycerol

The following examples consider a system composed a protein solvated by a mixture of water and glycerol, built with [Packmol](http://m3g.iqm.unicamp.br/packmol). The simulations were performed with [NAMD](https://www.ks.uiuc.edu/Research/namd/) with periodic boundary conditions and a NPT ensemble at room temperature and pressure. Molecular pictures were produced with [VMD](https://www.ks.uiuc.edu/Research/vmd/) and plots were produced with [Julia](https://julialang.org)'s [Plots](http://docs.juliaplots.org/latest/) library.

```@raw html
<center>
<img width=50% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/Data/system.png">
</center>
```
Image of the system of the example: a protein solvated by a mixture of glycreol (green) and water, at a concentration of 50%vv. 

### Index

- [Data](@ref data-example1)
- [MDDF, KB integrals, and group contributions](@ref)
- [2D density map](@ref)
- [3D density map](@ref)

## [Data](@id data-example1)

The [repository data](https://github.com/m3g/ComplexMixturesExamples/tree/main/Protein_in_Glycerol/Data) directory contains the a pdb file of the system (`system.pdb`) and a sample from the trajectory (`glyc50.dcd`), with a few frames. It also contains the result of running the `mddf` calculation on the complete trajectory, `results_glyc50.json`. :w

The sample trajectory is provided so that the first example can be run, yet do not expect that the results are the same, as the sampling is much lower in this case. The complete trajectory can be retrieved from [this link](https://drive.google.com/file/d/14M30jDHRwUM77hzbDphgbu8mcWFBcQrX/view?usp=sharing) (3GB file). 

## MDDF, KB integrals, and group contributions

Here we compute the minimum-distance distribution function, the Kirkwood-Buff integral, and the atomic contributions of the solvent to the density.
This example illustrates the regular usage of `ComplexMixtures`, to compute the minimum distance distribution function, KB-integrals and group contributions. 

```@raw html
<details><summary>Complete example code: click here!</summary>
```
```julia
# Activate environment in current directory
import Pkg; Pkg.activate(".")

# Run this once, to install necessary packages:
# Pkg.add(["ComplexMixtures", "PDBTools", "Plots", "LaTeXStrings"])

# Load packages
using ComplexMixtures
using PDBTools
using Plots, Plots.Measures
using LaTeXStrings

# The complete trajectory file can be downloaded from (3Gb):
# https://drive.google.com/file/d/14M30jDHRwUM77hzbDphgbu8mcWFBcQrX/view?usp=sharing

# The example output file is available at:
# 
# Load PDB file of the system
atoms = readPDB("./system.pdb")

# Select the protein and the GLYC molecules
protein = select(atoms, "protein")
glyc = select(atoms, "resname GLYC")

# Setup solute and solvent structures
solute = AtomSelection(protein, nmols=1)
solvent = AtomSelection(glyc, natomspermol=14)

# Path to the trajectory file
trajectory_file = "./glyc50_complete.dcd" 

# If the trajectory file is available, set run_mddf to "true" to 
# run the calculation from scratch
run_mddf = false
if run_mddf
    trajectory = Trajectory(trajectory_file, solute, solvent)
    results = mddf(trajectory)
    save(results, "glyc50_results.json")
else
    example_output = "./glyc50_results.json"
    results = load(example_output)
end

#
# Produce plots
#
# Default options for plots 
Plots.default(
    fontfamily="Computer Modern", 
    linewidth=2, 
    framestyle=:box, 
    label=nothing, 
    grid=false
)

#
# The complete MDDF and the Kirkwood-Buff Integral
#
plot(layout=(1, 2))
# plot mddf
plot!(results.d, results.mddf,
    xlabel=L"r/\AA", 
    ylabel="mddf", 
    subplot=1
)
hline!([1], linestyle=:dash, linecolor=:gray, subplot=1)
# plot KB integral
plot!(results.d, results.kb / 1000, #to L/mol
    xlabel=L"r/\AA", 
    ylabel=L"G_{us}/\mathrm{L~mol^{-1}}", 
    subplot=2
)
# size and margin
plot!(size=(800, 300), margin=4mm)
savefig("./mddf.png")

#
# Atomic contributions to the MDDF
#
hydroxyls = ["O1", "O2", "O3", "H1", "H2", "H3"]
aliphatic = ["C1", "C2", "HA", "HB", "HC", "HD"]
hydr_contrib = contributions(results, SolventGroup(hydroxyls))
aliph_contrib = contributions(results, SolventGroup(aliphatic))

plot(results.d, results.mddf, 
    xlabel=L"r/\AA", 
    ylabel="mddf", 
    size=(600, 400)
)
plot!(results.d, hydr_contrib, label="Hydroxyls")
plot!(results.d, aliph_contrib, label="Aliphatic chain")
hline!([1], linestyle=:dash, linecolor=:gray)
savefig("./mddf_atom_contrib.png")
```
```@raw html
</details><br>
```

#### Output 

The code above will produce the following plots, which contain the minimum-distance distribution of 
glycerol relative to the protein, and the corresponding KB integral:

```@raw html
<center>
<img width=100% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/MDDF/mddf.png">
</center>
```

and the same distribution function, decomposed into the contributions of the hydroxyl and aliphatic groups of glycerol:

```@raw html
<center>
<img width=70% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/MDDF/mddf_atom_contrib.png">
</center>
```

!!! note
    To change the options of the calculation, set the `Options` structure accordingly and pass it as a parameter to `mddf`. For example:
    ```julia
    options = Options(cutoff=10.)
    mddf(trajectory,options)
    ```
    The complete set of options available is described [here](@ref options).

## 2D density map

In this followup from the example aboave, we compute group contributions of the solute (the protein) to the MDDFs,
split into the contributions each protein residue. This allows the observation of the penetration of the solvent
on the structure, and the strength of the interaction of the solvent, or cossolvent, with each type of residue
in the structure.

```@raw html
<details><summary>Complete example code: click here!</summary>
```
```julia
# Activate environment in current directory
import Pkg; Pkg.activate(".")

# Run this once, to install necessary packages:
# Pkg.add(["ComplexMixtures", "PDBTools", "Plots", "LaTeXStrings"])

# Load packages
using ComplexMixtures
using PDBTools
using Plots, Plots.Measures
using LaTeXStrings

# Directory of the current script
script_dir = @__DIR__

# The complete trajectory file can be downloaded from (3Gb):
# https://drive.google.com/file/d/14M30jDHRwUM77hzbDphgbu8mcWFBcQrX/view?usp=sharing

# The example output file is available at:
# 
# Load PDB file of the system
atoms = readPDB("$script_dir/../Data/system.pdb")

# Select the protein and the GLYC molecules
protein = select(atoms, "protein")
glyc = select(atoms, "resname GLYC")

# Load example output file (computed in the previous script)
example_output = "./glyc50_results.json"
results = load(example_output)

#
# Plot a 2D map showing the contributions of some residues
#
residues = collect(eachresidue(protein))

# We will plot only the range 70:110, for clarity
irange = 70:110

# We create matrix of with a number of rows equal to the number
# of bins of the mddf histogram (length(results.d)) and a number of 
# columns equal to the number of residues
rescontrib = zeros(length(results.d), length(residues))

# Each column is then filled up with the contributions of each residue
for (ires, residue) in enumerate(residues)
    rescontrib[:, ires] .= contributions(results, SoluteGroup(residue))
end

# Plot only for distances within 1.5 and 3.5:
idmin = findfirst(d -> d > 1.5, results.d)
idmax = findfirst(d -> d > 3.5, results.d)

# Obtain pretty labels for the residues in the x-axis
xticks = PDBTools.residue_ticks(protein, first=70, last=110)

# Plot a contour courves with the density at each distance from
# each residue
contourf(irange, results.d[idmin:idmax], rescontrib[idmin:idmax, irange],
  color=cgrad(:tempo), linewidth=1, linecolor=:black,
  colorbar=:none, levels=5,
  xlabel="Residue", ylabel=L"r/\AA",
  xticks=xticks, xrotation=60,
  xtickfont=font(8, "Computer Modern"),
  size=(700, 400),
  margin=0.5Plots.PlotMeasures.cm
)
savefig("./density2D.png")
```
```@raw html
</details><br>
```

#### Output 

The code above will produce the following plot, which contains, for each residue, the contributions
of each residue to the distribution function of glycerol, within 1.5 to 3.5 $\mathrm{\AA}$ of the
surface of the protein.

```@raw html
<center>
<img width=70% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/Density2D/density2D.png">
</center>
```

## 3D density map

In this example we compute three-dimensional representations of the density map of Glycerol in the vicinity of a set of residues of a protein, from the minimum-distance distribution function. 

```@raw html
<details><summary>Complete example code: click here!</summary>
```
```julia
import Pkg; Pkg.activate(".")
using PDBTools
using ComplexMixtures

# PDB file of the system simulated
atoms = readPDB("./system.pdb")

# Load results of a ComplexMixtures run
results = load("./glyc50_results.json")

# Inform which is the solute
protein = select(atoms, "protein")
solute = AtomSelection(protein, nmols=1)

# Compute the 3D density grid and output it to the PDB file
# here we use dmax=3.5 such that the the output file is not too large
grid = grid3D(
    solute=solute,
    solute_atoms=protein,
    mddf_result=results,
    output_file="./grid.pdb",
    dmin=1.5,
    dmax=3.5
)
```
```@raw html
</details><br>
```

Here, the MDDF is decomposed at each distance according to the contributions of each *solute* (the protein) residue. The grid is created such that, at each point in space around the protein, it is possible to identify: 

1. Which atom is the closest atom of the solute to that point.

2. Which is the contribution of that atom (or residue) to the distribution function.

Therefore, by filtering the 3D density map at each distance one can visualize over the solute structure which are the regions that mostly interact with the solvent of choice at each distance. Typical images of such a density are:

```@raw html
<center>
<img width=100% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/Density3D/density3D_final.png">
</center>
```

In the figure on the left, the points in space around the protein are selected with the following properties: distance from the protein smaller than 2.0Å and relative contribution to the MDDF at the corresponding distance of at least 10% of the maximum contribution. Thus, we are selecting the regions of the protein corresponding to the most stable hydrogen-bonding interactions. The color of the points is the contribution to the MDDF, from blue to red. Thus, the most reddish-points corresponds to the regions where the most stable hydrogen bonds were formed. We have marked two regions here, on opposite sides of the protein, with arrows.

Clicking on those points we obtain which are the atoms of the protein contributing to the MDDF at that region. In particular, the arrow on the right points to the strongest red region, which corresponds to an Aspartic acid. These residues are shown explicitly under the density (represented as a transparent surface) on the figure in the center.

The figure on the right displays, overlapped with the hydrogen-bonding residues, the most important contributions to the second peak of the distribution, corresponding to distances from the protein between 2.0 and 3.5Å. Notably, the regions involved are different from the ones forming hydrogen bonds, indicating that non-specific interactions with the protein (and not a second solvation shell) are responsible for the second peak. 

### How to run this example:

Assuming that the input files are available in the script directory, just run the script with:

```bash
julia density3D.jl
```

Alternatively, open Julia and copy/paste or the commands in `density3D.jl` or use `include("./density3D.jl")`. These options will allow you to remain on the Julia section with access to the `grid` data structure that was generated and corresponds to the output `grid.pdb` file. 

This will create the `grid.pdb` file. [Here](https://raw.githubusercontent.com/m3g/ComplexMixturesExamples/main/Protein_in_Glycerol/Density3D/grid.vmd) we provide a previously setup VMD session that contains the data with the visualization choices used to generate the figure above. Load it with:

```bash
vmd -e grid.vmd
```

A short tutorial video showing how to open the input and output PDB files in VMD and produce images of the density is available here: 

```@raw html
<center>
<iframe width="560" style="height:315px" src="https://www.youtube.com/embed/V4Py44IKDh8" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</center>
```
