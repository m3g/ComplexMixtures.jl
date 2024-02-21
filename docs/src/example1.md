# Protein in water/glycerol

The following examples consider a system composed a protein solvated by a mixture of water and glycerol, built with [Packmol](http://m3g.iqm.unicamp.br/packmol). The simulations were performed with [NAMD](https://www.ks.uiuc.edu/Research/namd/) with periodic boundary conditions and a NPT ensemble at room temperature and pressure. Molecular pictures were produced with [VMD](https://www.ks.uiuc.edu/Research/vmd/) and plots were produced with [Julia](https://julialang.org)'s [Plots](http://docs.juliaplots.org/latest/) library.

```@raw html
<center>
<img width=50% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/Data/system.png">
</center>
```
Image of the system of the example: a protein solvated by a mixture of glycreol (green) and water, at a concentration of 50%vv. 

### Index

- [Data, packages, and execution](@ref data-example1)
- [MDDF, KB integrals, and group contributions](@ref mddf-example1)
- [2D density map](@ref 2D-map-example1)
- [3D density map](@ref 3D-map-example1)

## [Data, packages, and execution](@id data-example1)

The files required to run this example are:

- [system.pdb](https://raw.githubusercontent.com/m3g/ComplexMixturesExamples/main/Protein_in_Glycerol/Data/system.pdb): The PDB file of the complete system.
- [glyc50_traj.dcd](https://www.dropbox.com/scl/fi/zfq4o21dkttobg2pqd41m/glyc50_traj.dcd?rlkey=el3k6t0fx6w5yiqktyx96gzg6&dl=0): Trajectory file. This is a 1GB file, necessary for running from scratch the calculations.

To run the scripts, we suggest the following procedure:

1. Create a directory, for example `example1`.
2. Copy the required data files above to this directory.
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
   Alternativelly (and perhaps preferrably), copy line by line the content of the script into
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

## [2D density map](@id 2D-map-example1)

In this followup from the example aboave, we compute group contributions of the solute (the protein) to the MDDFs,
split into the contributions each protein residue. This allows the observation of the penetration of the solvent
on the structure, and the strength of the interaction of the solvent, or cossolvent, with each type of residue
in the structure.

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
<img width=70% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/Density2D/density2D.png">
</center>
```

## [3D density map](@id 3D-map-example1)

In this example we compute three-dimensional representations of the density map of Glycerol in the vicinity of a set of residues of a protein, from the minimum-distance distribution function. 

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
