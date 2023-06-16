# [Example](@id examples)

!!! note
    At [this repository](https://github.com/m3g/ComplexMixturesExamples) various examples are available illustrating the execution and possibilities of the package. Here we discuss one of these examples in detail.

The following examples consider a system composed a protein solvated by a mixture of water and glycerol, built with [Packmol](http://m3g.iqm.unicamp.br/packmol). The simulations were performed with [NAMD](https://www.ks.uiuc.edu/Research/namd/) with periodic boundary conditions and a NPT ensemble at room temperature and pressure. Molecular pictures were produced with [VMD](https://www.ks.uiuc.edu/Research/vmd/) and plots were produced with [Julia](https://julialang.org)'s [Plots](http://docs.juliaplots.org/latest/) library.

```@raw html
<center>
<img width=50% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/Data/system.png">
</center>
```

Image of the system of the example: a protein solvated by a mixture of glycreol (green) and water, at a concentration of 50%vv. 

## How to run this example

- Download and install [Julia](https://julialang.org)

- Install the required packages. Within Julia, do:
```julia-repl
julia> import Pkg

julia> Pkg.add(["ComplexMixtures", "PDBTools", "Plots", "LaTeXStrings", "Formatting"])
```

- Get the files:
```bash
git clone https://github.com/m3g/ComplexMixturesExamples
```

The files associated to the following examples are distributed at [this page](https://github.com/m3g/ComplexMixturesExamples/tree/main/Protein_in_Glycerol). 

## Data

The [Data](https://github.com/m3g/ComplexMixturesExamples/tree/main/Protein_in_Glycerol/Data) directory contains the a pdb file of the system (`system.pdb`) and a sample from the trajectory (`glyc50.dcd`), with a few frames. It also contains the result of running the `mddf` calculation on the complete trajectory, `results_glyc50.json`. This last file was produced by `ComplexMixtures`, as indicated in the following examples. 

The sample trajectory is provided so that the first example can be run, yet do not expect that the results are the same, as the sampling is much lower in this case. The complete trajectory can be retrieved from [this link](https://drive.google.com/file/d/14M30jDHRwUM77hzbDphgbu8mcWFBcQrX/view?usp=sharing) (3GB file). 

## Minimum-Distance Distribuion function

Here we compute the minimum-distance distribution function, the Kirkwood-Buff integral, and the atomic contributions of the solvent to the density.

This example illustrates the regular usage of `ComplexMixtures`, to compute the minimum distance distribution function, KB-integrals and group contributions. 

### How to run this example

```bash
cd ComplexMixturesExamples/Protein_in_Glycerol/MDDF
julia -t auto mddf.jl
```

### Detailed explanation of the example:

Loading the packages required for computing the MDDF.  

```julia
using PDBTools
using ComplexMixtures
```

Load the pdb file of the system using `PDBTools`:
```julia
atoms = readPDB("../Data/system.pdb")
```

Create arrays of atoms with the protein and Glycerol atoms, using the `select` function of the `PDBTools` package:
```julia
protein = select(atoms,"protein")
glyc = select(atoms,"resname GLYC")
```

Setup solute and solvent structures, required for computing the MDDF, with `Selection` function of the `ComplexMixtures` package:
```julia
solute = Selection(protein,nmols=1)
solvent = Selection(glyc,natomspermol=14)
```

Read and setup the Trajectory structure required for the computations:
```julia
trajectory = Trajectory("../Data/glyc50_complete.dcd",solute,solvent)
```

Run the calculation and get results:
```julia
results = mddf(trajectory)
```

!!! note
    To change the options of the calculation, set the `Options` structure accordingly and pass it as a parameter to `mddf`. For example:
    ```julia
    options = Options(cutoff=10.)
    mddf(trajectory,options)
    ```
    The complete set of options available is described [here](@ref options).


Save the reults to recover them later if required
```julia
save(results,"./glyc50.json")
```

The trajectory that was loaded was for a toy-example. The complete trajectory is available [here](https://drive.google.com/file/d/14M30jDHRwUM77hzbDphgbu8mcWFBcQrX/view?usp=sharing), but it is a 3GB file. The same procedure above was performed with that file and produced the `results_Glyc50.json` file, which is available in the Data directory here. We will continue with this file instead. 

Load the actual results obtained with the complete simulation:
```julia
results = load("../Data/results_glyc50.json")
```

Results are loaded, and now we can plot the data obtained.

### Produce plots

#### MDDF and Kirkwood-Buff integrals

Load some packages that we will use to produce the plots:
```julia
using Plots, Plots.PlotMeasures, LaTeXStrings
```

Some default options that make the plots prettier:
```julia
default(
    fontfamily="Computer Modern",
    linewidth=2, framestyle=:box, label=nothing, grid=false
)
```

First, we will plot the MDDF and the corresponding Kirkwood-Buff integral, which are available in the `results.mddf` and `results.kb` fields of the `results` data set. The distances are available in the `results.d` vector. We also plot here an horizontal line and save the figure as a `pdf` file.  

```julia
plot(layout=(1,2))
plot!(results.d,results.mddf,xlabel=L"r/\AA",ylabel="mddf",subplot=1)
hline!([1],linestyle=:dash,linecolor=:gray,subplot=1)
plot!(
    results.d,results.kb/1000, #to L/mol
    xlabel=L"r/\AA",ylabel=L"G_{us}/\mathrm{L~mol^{-1}}",
    subplot=2
)
plot!(size=(800,300),margin=4mm)
savefig("./mddf.pdf")
```

This will produce the following plot:

```@raw html
<center>
<img width=100% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/MDDF/mddf.png">
</center>
```

### Atomic contributions to the MDDF

Selecting the atoms corresponding to the hydroxyl groups, and of the aliphatic carbons of Glycerol. Here we list the types of the atoms as specified by the force-field.
```julia
hydroxyls = ["O1","O2","O3","H1","H2","H3"]
aliphatic = ["C1","C2","HA","HB","HC","HD"]
```

The `contributions` function of `ComplexMixtures` will extract from the result the contributions of each set of atoms to the total MDDF:

```julia
hydr_contributions = contributions(solvent,results.solvent_atom,hydroxyls)
aliph_contributions = contributions(solvent,results.solvent_atom,aliphatic)
```

And, finally, here we plot these group contributions on top of the total MDDF:

```julia
plot(results.d,results.mddf,xlabel=L"r/\AA",ylabel="mddf",size=(600,400))
plot!(results.d,hydr_contributions,label="Hydroxils")
plot!(results.d,aliph_contributions,label="Aliphatic chain")
hline!([1],linestyle=:dash,linecolor=:gray)
savefig("./mddf_atom_contrib.pdf")
```

This will produce the following figure:

```@raw html
<center>
<img width=70% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/MDDF/mddf_atom_contrib.png">
</center>
```

Note how hydroxyl clearly are the sole contribution to the peak at ~1.9 Angstroms, corresponding to hydrogen-bonding interactions. The aliphatic groups contribute importantly to the shoulder at larger distances, which correspond to non-specific interactions. 


## 2D residue contribution density map

In this example we compute the density map of Glycerol in the vicinity of a set of residues of a protein, from the minimum-distance distribution function. 

The MDDF can be decomposed in the contributions of each atom of the solute or of the solvent. Here, we sum up te contributions of all the atoms of each residue of the solute, which is a protein, and plot a density map with the final information. The output figure obtained is:

```@raw html
<center>
<img width=70% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/Density2D/density2D.png">
</center>
```

### How to run this example:

```bash
cd ComplexMixturesExamples/Protein_in_Glycerol/Density2D
julia density2D.jl
```

### Detailed explanation of the example:

Here, we use the `contourf` function of the `Plots` package of Julia. A detailed explanation of the input file `density2D.jl` is provide below: 

### Loading packages that will be used:


```julia
using Plots
using LaTeXStrings
using Formatting
using ComplexMixtures, PDBTools
```

### Some default options so the plot looks nice
```julia
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2, framestyle=:box, label=nothing
)
```

### Read the PDB file (using PDBTools)
```julia
pdb = readPDB("./system.pdb")
```

### Load results of the ComplexMixtures run
```julia
R = load("./results_glyc50.json")  
```

### Define which are the solute molecules (the protein)
```julia
protein = select(pdb,"protein")
solute = Selection(protein,nmols=1)
```

### Define which are the solvent molecules (Glycerol here)
```julia
glycerol = select(pdb,"resname GLYC")
solvent = Selection(glycerol,natomspermol=14)
```

### Retrive the resiude contribution data

Collect which are the protein residues 
```julia
residues = collect(eachresidue(protein))
```

Set a matrix that will store the results, with a number of lines corresponding to the length of the MDDF histogram, and with a number of columns corresponding to the number of residues:
```julia
rescontrib = zeros(length(R.mddf),length(residues))
```

Now, collect the contribution of each residue as a column of the above matrix. The notation `pairs(residues)` returns tuples containg the index `ires` and the corresponding residue. The `.=` symbol sets each element of the corresponding column of the  `rescontrib` matrix to the output of `contributions` (by broadcasting).  
```julia
for (ires,residue) in pairs(residues)
  rescontrib[:,ires] .= contributions(solute,R.solute_atom,residue)
end
```

### Plot only for distances within 1.5 and 3.5:

Here, we will plot only the contributions from residue `70` to residue `110`, and from distances ranging from `1.5` to `3.5` which is where most of the action occurs:
```julia
irange=70:110
idmin = findfirst( d -> d > 1.5, R.d)
idmax = findfirst( d -> d > 3.5, R.d)
```

To obtain pretty labels for the residues in the x-axis, we retrieve the one-letter residue names and concatenate them with the residue number converted to strings:

```julia
labels = PDBTools.oneletter.(resname.(residues)).*format.(resnum.(residues))
```

And, finally, we produce the plot, with a series of options that make this particular contour plot look nice:

```julia
contourf(
    irange, # x
    R.d[idmin:idmax], # y
    rescontrib[idmin:idmax,irange], # z
    xlabel="Residue", ylabel=L"r/\AA",
    xticks=(irange,labels[irange]), xrotation=60,
    xtickfont=font(6,plot_font),
    color=cgrad(:tempo), linewidth=0.1, linecolor=:black,
    colorbar=:none, levels=5,
    size=(500,280)
)
```

The final figure is saved as a `pdf` file:
```julia
savefig("./density2D.pdf")
```

# 3D residue contribution density map

In this example we compute three-dimensional representations of the density map of Glycerol in the vicinity of a set of residues of a protein, from the minimum-distance distribution function. 

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

A short tutorial video showing how to open the input and output PDB files in VMD and produce images of the density is available here: 

```@raw html
<center>
<iframe width="560" style="height:315px" src="https://www.youtube.com/embed/V4Py44IKDh8" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
</center>
```

### How to run this example:

```bash
cd ComplexMixturesExamples/Protein_in_Glycerol/Density3D
julia density3D.jl
```

Alternatively, open Julia and copy/paste or the commands in `density3D.jl` or use `include("./density3D.jl")`. These options will allow you to remain on the Julia section with access to the `grid` data structure that was generated and corresponds to the output `grid.pdb` file. 

This will create (actually overwrite) the `grid.pdb` file. Here we provide a previously setup VMD session that contains the data with the visualization choices used to generate the figure above. Load it with:

```bash
vmd -e grid.vmd
``` 

### Detailed explanation of the example:

Initially we load the `ComplexMixtures` and `PDBTools` packages:

```julia
using ComplexMixtures, PDBTools
```

With the `readPDB` function of `PDBTools`, we read the  PDB file of the system simulated:
```julia
pdb = readPDB("../Data/system.pdb")
```

and using `ComplexMixtures`, we load the results from the calculation of the MDDF of Glycerol around the protein, which was computed previously:
```julia
R = load("../Data/results_glyc50.json")  
```

The solute here is the protein, and we need to setup the structures that define which atoms and type of solute it is. First, we select from the atoms of the pdb file of the system, those belonging to the protein, using `select` from `PDBTools`:
```julia
protein = select(pdb,"protein")
```

and then we define the solute structure that is actually used in `ComplexMixtures`, by passing those atoms and specifying that the solute is a single molecule to the `Selection` function of `ComplexMixtures`:
```julia
solute = Selection(protein,nmols=1)
```

The 3D grid representing the density around the protein is computed with the `grid3D` function provided by `ComplexMixtures`. It receives the `solute` structure (of type `Selection`), the list of solute atoms (of type `PDBTools.Atoms`, as the `protein` selection above), the name of the output file and some optional parameters to define the grid. Here we compute the grid only between 1.5 and 3.5Å, characterizing the first and second solvation shells. The grid has by default a `step` of 0.5Å. 

```julia
grid = grid3D(
    solute=solute,
    solute_atoms=protein,
    mddf_result=R,
    output_file="grid.pdb",
    dmin=1.5,
    dmax=3.5
)
```

The command above will generate the grid, save it to `grid.pdb` and let it available in the `grid.pdb` array of atoms, for further inspection, if desired. 

By changing `dmin`, `dmax`, and `step`, one controls the grid size and resolution. This may generate very large output files.
                                                                                                      

































