# [Protein in water/glycerol](@id example1)

The following examples consider a system composed a protein solvated by a mixture of water and glycerol, built with [Packmol](http://m3g.iqm.unicamp.br/packmol). The simulations were performed with [NAMD](https://www.ks.uiuc.edu/Research/namd/) with periodic boundary conditions and a NPT ensemble at room temperature and pressure. Molecular pictures were produced with [VMD](https://www.ks.uiuc.edu/Research/vmd/) and plots were produced with [Julia](https://julialang.org)'s [Plots](http://docs.juliaplots.org/latest/) library.

```@raw html
<center>
<img width=50% src="../figures/prot_glyc_system.png">
</center>
```
Image of the system of the example: a protein solvated by a mixture of glycerol (green) and water, at a concentration of 50%vv. 

### Index

- [Data, packages, and execution](@ref data-example1)
- [MDDF, KB integrals, and group contributions](@ref mddf-example1)
- [2D density map](@ref 2D-map-example1)
- [3D density map](@ref 3D-map-example1)

## [Data, packages, and execution](@id data-example1)

The files required to run this example are available at [[this link]](https://zenodo.org/records/17202344/files/example1.zip?download=1), and are: 

- `system.pdb`: The PDB file of the complete system.
- `glyc50_traj.dcd`: Trajectory file. This is a 1GB file, necessary for running from scratch the calculations.

To run the scripts, we suggest the following procedure:

1. Create a directory, for example `example1`.
2. Unzip and copy the required data files above to this directory.
3. Launch `julia` in that directory, activate the directory environment, and install the required packages. 
   This is done by launching Julia and executing:
   ```julia
   import Pkg 
   Pkg.activate(".")
   Pkg.add(["ComplexMixtures", "PDBTools", "Plots", "LaTeXStrings", "EasyFit"])
   exit()
   ```
4. Copy the code of each script in to a file, and execute with:
   ```julia
   julia -t auto script.jl
   ```
   Alternatively (and perhaps preferably), copy line by line the content of the script into
   the Julia REPL, to follow each step of the calculation. For a more advanced Julia usage,
   we suggest the [VSCode IDE](https://code.visualstudio.com/) with the 
   [Julia Language Support](https://www.julia-vscode.org/docs/dev/gettingstarted/) extension. 

## [MDDF, KB integrals, and group contributions](@id mddf-example1)

Here we compute the minimum-distance distribution function, the Kirkwood-Buff integral, and the atomic contributions of the solvent to the density.
This example illustrates the regular usage of `ComplexMixtures`, to compute the minimum distance distribution function, KB-integrals and group contributions. 

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example1/script1.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

#### Output 

The code above will produce the following plots, which contain the minimum-distance distribution of 
glycerol relative to the protein, and the corresponding KB integral:

```@raw html
<center>
<img width=100% src="../assets/scripts/example1/mddf.png">
</center>
```

and the same distribution function, decomposed into the contributions of the hydroxyl and aliphatic groups of glycerol:

```@raw html
<center>
<img width=70% src="../assets/scripts/example1/mddf_atom_contrib.png">
</center>
```

!!! note
    To change the options of the calculation, set the `Options` structure accordingly and pass it as a parameter to `mddf`. For example:
    ```julia
    options = Options(bulk_range=(10.0, 15.0), stride=5)
    mddf(trajectory_file, solute, solvent, options)
    ```
    The complete set of options available is described [here](@ref options).

## [2D density map](@id 2D-map-example1)

In this followup from the example above, we compute group contributions of the solute (the protein) to the MDDFs,
split into the contributions each protein residue. This allows the observation of the penetration of the solvent
on the structure, and the strength of the interaction of the solvent, or cosolvent, with each type of residue
in the structure. The `ResidueContributions` and `Plots.contourf` auxiliary functions, [documented here](@ref 2D_per_residue), are used:  

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example1/script2.jl", String))
\`\`\`
""")
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
<img width=70% src="../assets/scripts/example1/density2D.png">
</center>
```

## [3D density map](@id 3D-map-example1)

In this example we compute three-dimensional representations of the density map of Glycerol in the vicinity of a set of residues of a protein, from the minimum-distance distribution function. 
For further information about the `grid3D` function, go to the [3D grid map](@ref grid3D) section of the manual.

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example1/script3.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

Here, the MDDF is decomposed at each distance according to the contributions of each *solute* (the protein) atom. The grid is created such that, at each point in space around the protein, it is possible to identify: 

1. Which atom is the closest atom of the solute to that point.

2. Which is the contribution of that atom (or residue) to the distribution function.

Therefore, by filtering the 3D density map at each distance one can visualize over the solute structure which are the regions that mostly interact with the solvent of choice at each distance. The views are created with the `visualize` function of [PDBTools](https://m3g.github.io/PDBTools.jl/stable/visualization/) (version 3.41.0 or greater), and the script saves them as the `grid_hbonds.html` and `grid_second_peak.html` files, which can be opened in any web browser. Typical views are shown below (rotate, zoom, and hover over the points to identify the protein atoms):

```@example example1-grid3D
using ComplexMixtures, PDBTools # hide
using ComplexMixtures: data_dir # hide
atoms = read_pdb(joinpath(data_dir, "NAMD", "Protein_in_Glycerol", "system.pdb")) # hide
protein = select(atoms, "protein") # hide
results = load(joinpath(data_dir, "NAMD", "Protein_in_Glycerol", "protein_glyc.json")) # hide
grid = grid3D(results, atoms; dmin=1.5, dmax=3.5, silent=true) # hide
hbonds = filter(p -> occup(p) < 2.0 && beta(p) > 0.1, grid) # hide
visualize( # hide
    protein => (color="white",), # hide
    hbonds => (style=:dots, color_by=beta.(hbonds), colormap=:rwb, color_range=(0.4, -0.4)), # hide
    protein => (selection="resname GLU ASP", style=:ballandstick, color="green"); # hide
    height=500, # hide
) # hide
```

In the view above, the points in space around the protein are selected with the following properties: distance from the protein smaller than 2.0Å and relative contribution to the MDDF at the corresponding distance of at least 10% of the maximum contribution. Thus, we are selecting the regions of the protein corresponding to the most stable hydrogen-bonding interactions. The protein is shown as a white cartoon, the Glutamic and Aspartic acid residues as green balls and sticks, and the grid points as dots, colored by the contribution to the MDDF, from white to red. Thus, the most reddish points correspond to the regions where the most stable hydrogen bonds were formed. 

Hovering over those points we obtain which are the atoms of the protein contributing to the MDDF at that region (for example, `ASP14:P OD1` for the `OD1` atom of residue 14). In particular, the strongest red regions correspond to the carboxylate groups of Aspartic and Glutamic acid residues. 

The view below displays the most important contributions to the second peak of the distribution, corresponding to distances from the protein between 2.0 and 3.5Å, with the Threonine, Asparagine, and Tryptophan residues shown as balls and sticks:

```@example example1-grid3D
second_peak = filter(p -> 2.0 < occup(p) < 3.5 && beta(p) > 0.1, grid) # hide
visualize( # hide
    protein => (color="white",), # hide
    second_peak => (style=:dots, color_by=beta.(second_peak), colormap=:rwb, color_range=(0.4, -0.4)), # hide
    protein => (selection="resname THR ASN TRP", style=:ballandstick); # hide
    height=500, # hide
) # hide
```

Notably, the regions involved are different from the ones forming hydrogen bonds, indicating that non-specific interactions with the protein (and not a second solvation shell) are responsible for the second peak. 

!!! tip
    Other ways to project the contributions into the structure, such as coloring the atoms or residues of the protein 
    by their contributions to coordination numbers or KB integrals, are shown [here](@ref grid3D-visualization).

### How to run this example:

Assuming that the input files are available in the script directory, just run the script with:

```bash
julia -t auto script3.jl
```

Alternatively, open Julia and copy/paste the commands of `script3.jl` or use `include("./script3.jl")`. These options will allow you to remain on the Julia session with access to the `grid` data structure that was generated, and to display the views interactively (in VSCode or in the browser, if running from the REPL).

### Visualization with VMD

The grid can also be written to a PDB file and visualized with other software. An example of the visualization of the grid of this example with [VMD](https://www.ks.uiuc.edu/Research/vmd/) is available [here](@ref grid3D-vmd).

### References

The interactive 3D views were rendered with [3Dmol.js](https://3dmol.csb.pitt.edu), through the `visualize` function of [PDBTools](https://m3g.github.io/PDBTools.jl/stable/visualization/): N. Rego and D. Koes, 3Dmol.js: molecular visualization with WebGL, Bioinformatics 31, 1322-1324 (2015). [https://doi.org/10.1093/bioinformatics/btu829](https://doi.org/10.1093/bioinformatics/btu829)
