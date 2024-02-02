# Polyacrylamide in DMDF

In this example we illustrate how the solvation structure of a polymer can be studied with ComplexMixtures.jl. The system is a 5-mer segment of polyacrylamide (PAE - capped with methyl groups), solvated with dimethylformamide (DMF). The system is interesting because of the different functional groups and polarities involved in the interactions of DMF with PAE. A snapshot of the system is shown below.

```@raw html
<center>
<img width=30% src="https://raw.githubusercontent.com/m3g/ComplexMixturesExamples/main/Polyacrylamide_in_DMF/system.png">
</center>
```

The structures of DMF and of the polyacrylamide segment are:

```@raw html
<center>
<table><tr>
<td><img src="https://raw.githubusercontent.com/m3g/ComplexMixturesExamples/main/Polyacrylamide_in_DMF/simulation/dmf.png" height=100px></td>
<td><img src="https://raw.githubusercontent.com/m3g/ComplexMixturesExamples/main/Polyacrylamide_in_DMF/simulation/polyacrylamide.png" height=100px></td>
</tr>
<tr>
<td align=center>DMF</td>
<td align=center>Polyacrylamide</td>
</tr>
</table>
</center>
```

### Index

- [Data, packages, and execution](@ref data-example2)
- [MDDF and KB integrals](@ref mddf-example2)
- [Group contributions](@ref groups-example2)
- [2D density map](@ref 2Dmap-example2)

## [Data and execution](@id data-example2)

The files required to run this example are:

- [equilibrated.pdb](https://raw.githubusercontent.com/m3g/ComplexMixturesExamples/main/Polyacrylamide_in_DMF/simulation/equilibrated.pdb): The PDB file of the complete system.
- [traj_Polyacry.dcd](https://www.dropbox.com/scl/fi/jwafhgxaxuzsybw3y8txd/traj_Polyacry.dcd?rlkey=p4bn65m0pkuebpfm0hf158cdm&dl=0): Trajectory file. This is a 275Mb file, necessary for running from scratch the calculations.

To run the scripts, we suggest the following procedure:

1. Create a directory, for example `example2`.
2. Copy the required data files above to this directory.
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
   Alternativelly (and perhaps preferrably), copy line by line the content of the script into
   the Julia REPL, to follow each step of the calculation.

## [MDDF and KB integrals](@id mddf-example2)

Here we compute the minimum-distance distribution function, the Kirkwood-Buff integral, and the atomic contributions of the solvent to the density.
This example illustrates the regular usage of `ComplexMixtures`, to compute the minimum distance distribution function, KB-integrals and group contributions. 

```@raw html
<details><summary>Complete example code: click here!</summary>
```
```julia
import Pkg; Pkg.activate(".")

using PDBTools
using ComplexMixtures
using Plots
using LaTeXStrings
using EasyFit: moveavg

# The full trajectory file is available at: 
# https://www.dropbox.com/scl/fi/jwafhgxaxuzsybw3y8txd/traj_Polyacry.dcd?rlkey=p4bn65m0pkuebpfm0hf158cdm&dl=0 
trajectory_file = "./traj_Polyacry.dcd"

# Load a PDB file of the system
system = readPDB("./equilibrated.pdb")

# Select the atoms corresponding DMF molecules
dmf = select(system, "resname DMF")

# Select the atoms corresponding to the Poly-acrylamide
acr = select(system, "resname FACR or resname ACR or resname LACR")

# Set the solute and the solvent selections for ComplexMixtures
solute = AtomSelection(acr, nmols=1)
solvent = AtomSelection(dmf, natomspermol=12)

# Set the trajectory structure
trajectory = Trajectory(trajectory_file, solute, solvent)

# Use a large dbulk distance for better KB convergence
options = Options(dbulk=20.)

# Compute the mddf and associated properties
results = mddf(trajectory, options)

# Save results to file for later use
save(results, "./mddf.json")
println("Results saved to ./mddf.json file")

# Plot the MDDF and KB integrals
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=false,
    palette=:tab10
)
scalefontsizes(); scalefontsizes(1.3)

# Plot the MDDF of DMF relative to PolyACR and its corresponding KB integral
mddf_smoothed = movavg(results.mddf,n=10).x # Smooth example with a running average
plot(layout=(2,1))
plot!(results.d, mddf_smoothed, ylabel="MDDF", subplot=1, xlims=(0,20))

smoothed_kb = movavg(results.kb,n=10).x # smooth kb
plot!(results.d, smoothed_kb, 
    xlabel=L"\textrm{Distance / \AA}",
    ylabel=L"\textrm{KB~/~cm^3~mol^{-1}}",
    xlim=(0,20),
    subplot=2
)
savefig("./mddf_kb.png")
println("Plot saved to mddf_kb.png")
```
```@raw html
</details><br>
```

#### Output 

The distribution of DMF molecules around polyacrylamide is shown below. There is a peak at ~2.5Angs, indicating favorable non-specific interactions between the solvent molecules and the polymer. The peak is followed by a dip and diffuse peaks at higher distances. Thus, the DMF molecules are structured around the polymer, but essentially only in the first solvation shell.  

![](https://raw.githubusercontent.com/m3g/ComplexMixturesExamples/main/Polyacrylamide_in_DMF/results/mddf_kb.png)

The KB integral in a bicomponent mixture converges to the (negative of the) apparent molar volume of the solute. It is negative, indicating that the accumulation of DMF in the first solvation shell of the polymer is not enough to compensate the excluded volume of the solute. 

## [Group contributions](@id groups-example2)

voltar

#### Output 

## [2D density map](@id 2Dmap-example2)

voltar

#### Output 

The code above will produce the following plot, which contains, for each residue, the contributions
of each residue to the distribution function of glycerol, within 1.5 to 3.5 $\mathrm{\AA}$ of the
surface of the protein.

```@raw html
<center>
<img width=70% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/Density2D/density2D.png">
</center>
```