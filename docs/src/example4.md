# Glycerol/water mixture

This example illustrates the use of ComplexMixtures.jl to study the solution structure of a crowded (1:1 molar fraction) solution of glycerol in water. Here, we compute the distribution function and atomic contributions associated to the inter-species interactions (water-glycerol) and the glycerol-glycerol auto-correlation function. This example aims to illustrate how to obtain a detailed molecular picture of the solvation structures in an homogeneous mixture.

The system simulated consists of 1000 water molecules (red) and 1000 glycerol molecules (purple).
```@raw html
<center>
<img src="../figures/glyc_wat_system.png" width=30%>
</center>
```

### Index

- [Data, packages, and execution](@ref data-example4)
- [Glycerol-Glycerol and Water-Glycerol distribution functions](@ref glyc_mddf-example4)
- [Glycerol group contributions to MDDFs](@ref glyc-groups-example4)
- [2D map of group contributions](@ref map-example4)
- [3D density map of water around glycerol](@ref 3Dmap-example4)

## [Data, packages, and execution](@id data-example4)

The files required to run this example are available at [[this link]](https://zenodo.org/records/17202344/files/example4.zip?download=1), and are: 

- `equilibrated.pdb`: The PDB file of the complete system.
- `traj_Glyc.dcd`: Trajectory file. This is a 200Mb file, necessary for running from scratch the calculations.

To run the scripts, we suggest the following procedure:

1. Create a directory, for example `example4`.
2. Unzip and copy the required data files above to this directory.
3. Launch `julia` in that directory: activate the directory environment, and install the required packages. This launching Julia and executing:
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
   the Julia REPL, to follow each step of the calculation.

## [Glycerol-Glycerol and Water-Glycerol distribution functions](@id glyc_mddf-example4)

The first and most simple analysis is the computation of the minimum-distance distribution functions between the components of the solution. In this example we focus on the distributions of the two components relative to the glycerol molecules. Thus, we display the glycerol auto-correlation function, and the water-glycerol correlation function in the first panel of the figure below. The second panel displays the KB integrals of the two components computed from each of these distributions.

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example4/script1.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

![](./assets/scripts/example4/mddf_kb.png)

Both water and glycerol form hydrogen bonds with (other) glycerol molecules, as indicated by the peaks at ~1.8$$\mathrm{\AA}$$. The auto-correlation function of glycerol shows a more marked second peak corresponding to non-specific interactions, which (as we will show) are likely associated to interactions of its aliphatic groups.

The KB integrals in the second panel show similar values for water and glycerol, with the KB integral for water being slightly greater. This means that glycerol molecules are (slightly, if the result is considered reliable) preferentially hydrated from a macroscopic standpoint.

## [Glycerol group contributions to MDDFs](@id glyc-groups-example4)

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example4/script2.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

![](./assets/scripts/example4/mddf_group_contributions.png)


## [2D map of group contributions](@id map-example4)

The above distributions can be split into the contributions of each glycerol chemical group. The 2D maps below display this decomposition.

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example4/script3.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

![](./assets/scripts/example4/GlycerolWater_map.png)

The interesting result here is that the $$\mathrm{CH}$$ group of glycerol is protected from both solvents. There is a strong density augmentation at the vicinity of hydroxyl groups, and the second peak of the MDDFs is clearly associated to interactions with the $$\mathrm{CH_2}$$ groups.

## [3D density map of water around glycerol](@id 3Dmap-example4)

The contributions of the glycerol atoms to the MDDF of water can also be projected into the space around 
a glycerol molecule, with the [`grid3D`](@ref grid3D) function, and visualized with the `visualize` function of 
[PDBTools](https://m3g.github.io/PDBTools.jl/stable/visualization/) (version 3.41.0 or greater). Since the solute
(glycerol) is composed of many molecules, the grid is built around one of them, and the contributions of each
atom are averaged over all glycerol molecules.

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example4/script4.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

#### Output 

The script saves the views as the `grid_hbonds.html` and `isosurfaces.html` files, which can be opened in any 
web browser (rotate, zoom, and hover over the points to identify the glycerol atoms closest to each point). 

The first view displays the regions associated with the hydrogen-bonding peak of the MDDF (distances 
smaller than 2.0Å), with at least 10% of the maximum contribution. The glycerol molecule is shown in balls 
and sticks, and the grid points as dots, colored from white to red according to their contribution to the MDDF: 

```@raw html
<center>
<iframe src="../assets/scripts/example4/grid_hbonds.html" style="width: 100%; height: 450px; border: none;"></iframe>
</center>
```

In the second view, the grid is converted into volumetric data with the [`volumetric_data`](@ref) function, and 
smoothed with a Gaussian function of width `sigma=0.5` Å. The isosurfaces correspond to 50% (orange, transparent) 
and 75% (red) of the maximum value of the smoothed data:

```@raw html
<center>
<iframe src="../assets/scripts/example4/isosurfaces.html" style="width: 100%; height: 450px; border: none;"></iframe>
</center>
```

The highest densities of water (red) are found at hydrogen-bonding distances from the hydroxyl groups of glycerol.
The volumetric data is also written to the `density.dx` file (OpenDX format), which can be visualized in other
software, such as [VMD](https://www.ks.uiuc.edu/Research/vmd/).

### References

The interactive 3D views were rendered with [3Dmol.js](https://3dmol.csb.pitt.edu), through the `visualize` function of [PDBTools](https://m3g.github.io/PDBTools.jl/stable/visualization/): N. Rego and D. Koes, 3Dmol.js: molecular visualization with WebGL, Bioinformatics 31, 1322-1324 (2015). [https://doi.org/10.1093/bioinformatics/btu829](https://doi.org/10.1093/bioinformatics/btu829)
