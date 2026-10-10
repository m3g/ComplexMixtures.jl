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
- [Preferential interactions and m-values](@ref preferential-example1)
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
    options = Options(cutoff=15.0, stride=5)
    mddf(trajectory_file, solute, solvent, options)
    ```
    The complete set of options available is described [here](@ref options).

## [Preferential interactions and m-values](@id preferential-example1)

The KB integral of a single component of the solvent is difficult to interpret by itself: for a large 
solute, it is dominated by the excluded volume of the solute, and is thus large and negative for all
components (about -20 L mol⁻¹ for both water and glycerol here). The thermodynamically relevant quantity
is the *difference* between the KB integrals of the cosolvent and of water, which defines the 
**preferential interaction parameter** of the cosolvent (``c``, glycerol) relative to water (``w``), 
around the protein (``p``) [1,2]:

```math
\Gamma_{pc} = \rho_c \left( G_{pc} - G_{pw} \right)
```

where ``\rho_c`` is the bulk molar concentration of the cosolvent. ``\Gamma_{pc} > 0`` indicates that the 
cosolvent accumulates on the protein surface relative to water (preferential solvation by the cosolvent),
and ``\Gamma_{pc} < 0`` that it is preferentially excluded. The excluded volume contributions cancel in the
difference, and ``\Gamma_{pc}`` can be interpreted as the excess number of cosolvent molecules
in the solvation domain of the protein, relative to the composition of the bulk solution.

``\Gamma_{pc}`` is related to the dependence of the chemical potential of the protein, ``\mu_p``, on the 
activity of the cosolvent, ``a_c``, at constant temperature and pressure:

```math
\left(\frac{\partial \mu_p}{\partial \ln a_c}\right)_{T,P} = -RT\,\Gamma_{pc}
```

and, therefore, on the concentration of the cosolvent [1]:

```math
\left(\frac{\partial \mu_p}{\partial \rho_c}\right)_{T,P} = 
-\frac{RT \left( G_{pc} - G_{pw} \right)}{1 + \rho_c\left( G_{cc} - G_{cw} \right)}
```

The denominator, ``\left(\partial \ln a_c / \partial \ln \rho_c\right)^{-1}``, accounts for the 
non-ideality of the water-cosolvent mixture, and depends on the KB integrals between the 
components of the solvent, ``G_{cc}`` and ``G_{cw}``. If the mixture is ideal, ``G_{cc} = G_{cw}``,
and the denominator is one. 

The **m-value** of a process involving the protein (for example, unfolding) is the derivative of its
free energy with respect to the cosolvent concentration, and is thus obtained from the difference of
the quantities above in the two states (for example, unfolded, U, and folded, F, states):

```math
m = \frac{\partial \Delta G_{F \to U}}{\partial \rho_c} = 
-\frac{RT \left[ \left( G_{pc} - G_{pw} \right)^U - \left( G_{pc} - G_{pw} \right)^F \right]}{1 + \rho_c\left( G_{cc} - G_{cw} \right)}
```

The simulation of a single state, as here, provides the solvation contribution of that state
to the m-value.

### Computing the preferential interaction parameter

To compute ``\Gamma_{pc}``, the KB integral of water relative to the protein is required, in addition 
to that of glycerol, computed above. The following script computes the MDDF of water from the same
trajectory, and computes ``\Gamma_{pc}`` and the derivative of the chemical potential of the protein, 
assuming that the water-glycerol mixture is ideal:

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example1/script4.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

#### Output

The KB integrals of glycerol and water relative to the protein, and the corresponding preferential 
interaction parameter, as a function of the distance to the protein, are:

```@raw html
<center>
<img width=100% src="../assets/scripts/example1/preferential_interaction.png">
</center>
```

The KB integrals, and thus ``\Gamma_{pc}``, are stable for distances greater than about 8 Å, indicating that
the bulk solution is reached at these distances from the protein surface. The script prints the values
averaged between 10 and 15 Å:

```
Bulk concentration of glycerol: 5.95 mol/L
G_pc - G_pw = 2.1 L/mol
Preferential interaction parameter, Γ = 12.6
∂μ_p/∂ln(a_c) = -RTΓ = -7.4 kcal/mol
∂μ_p/∂ρ_c (ideal mixture) = -1.3 kcal mol⁻¹ M⁻¹
```

In this simulation, thus, glycerol is preferentially bound to the protein: there are about 13 more glycerol 
molecules in the solvation domain of the protein than expected from the bulk composition. The chemical potential
of the protein decreases with the increase of the glycerol concentration (assuming the ideality of the
water-glycerol mixture). The negative values of ``\Gamma_{pc}`` at short distances reflect the fact that water,
being smaller, approaches the protein surface more closely than glycerol. The 2D and 3D density maps below
illustrate which regions of the protein are responsible for the interactions with glycerol.

!!! note
    The sign and magnitude of ``\Gamma_{pc}`` depend on the force field and on the sampling. Here, a 
    single trajectory of a relatively small system was used, and the results are illustrative. 

### Accounting for the non-ideality of the solvent

The non-ideality factor, ``1 + \rho_c(G_{cc} - G_{cw})``, can be obtained from experimental 
activity data of the water-cosolvent mixture, or from a simulation of the water-cosolvent mixture 
without the protein, at the same composition. In this case, the solvent components are small molecules, 
and each can be represented by a single atom (for example, the central carbon of glycerol, `C2`, 
and the oxygen atom of water, `OH2`). Thus, the KB integrals are computed from radial distribution functions,
and are obtained with the [`kbi`](@ref) function. For instance, assuming that 
`mixture.pdb` and `mixture.dcd` are the files of a simulation of the water-glycerol mixture:

```julia
using ComplexMixtures, PDBTools
atoms = read_pdb("mixture.pdb")
glyc_C2 = AtomSelection(select(atoms, "resname GLYC and name C2"), natomspermol=1)
water_O = AtomSelection(select(atoms, "water and name OH2"), natomspermol=1)
# glycerol-glycerol: solute and solvent are the same
R_cc = mddf("mixture.dcd", glyc_C2, Options(cutoff=20.0))
# glycerol-water 
R_cw = mddf("mixture.dcd", glyc_C2, water_O, Options(cutoff=20.0))
# KB integrals (L/mol), and the bulk concentration of glycerol (mol/L)
G_cc = kbi(R_cc) / 1000
G_cw = kbi(R_cw) / 1000
ρ_c = overview(R_cc).density.solvent_bulk
# Non-ideality factor, as a function of L (take the value where it is stable)
nonideality = @. 1 + ρ_c * (G_cc - G_cw)
```

The value of ``\partial \mu_p / \partial \rho_c`` computed above, assuming ideality, must then be divided by
this factor.

### References

1. S. Shimizu, Estimating hydration changes upon biomolecular reactions from osmotic stress, high pressure,
   and preferential hydration experiments. *Proc. Natl. Acad. Sci. USA* 101, 1195–1199 (2004).
   [DOI: 10.1073/pnas.0305836101](https://doi.org/10.1073/pnas.0305836101)
2. L. Martínez, S. Shimizu, Molecular interpretation of preferential interactions in protein solvation: 
   a solvent-shell perspective by means of minimum-distance distribution functions.
   *J. Chem. Theory Comput.* 13, 6358–6372 (2017). 
   [DOI: 10.1021/acs.jctc.7b00599](https://doi.org/10.1021/acs.jctc.7b00599)

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
