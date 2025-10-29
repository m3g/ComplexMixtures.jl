# Polyacrylamide in DMDF

In this example we illustrate how the solvation structure of a polymer can be studied with ComplexMixtures.jl. The system is a 5-mer segment of polyacrylamide (PAE - capped with methyl groups), solvated with dimethylformamide (DMF). The system is interesting because of the different functional groups and polarities involved in the interactions of DMF with PAE. A snapshot of the system is shown below.

```@raw html
<center>
<img width=30% src="../figures/poly_dmf_system.png">
</center>
```

The structures of DMF and of the polyacrylamide segment are:

```@raw html
<center>
<table><tr>
<td><img src="../figures/dmf.png" height=100px></td>
<td><img src="../figures/polyacrylamide.png" height=100px></td>
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

## [Data, packages, and execution](@id data-example2)

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
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example2/script1.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

#### Output 

The distribution of DMF molecules around polyacrylamide is shown below. There is a peak at ~2.5Angs, indicating favorable non-specific interactions between the solvent molecules and the polymer. The peak is followed by a dip and diffuse peaks at higher distances. Thus, the DMF molecules are structured around the polymer, but essentially only in the first solvation shell.  

![](./assets/scripts/example2/mddf_kb.png)

The KB integral in a bicomponent mixture converges to the (negative of the) apparent molar volume of the solute. It is negative, indicating that the accumulation of DMF in the first solvation shell of the polymer is not enough to compensate the excluded volume of the solute. 

## [Group contributions](@id groups-example2)

The MDDF can be decomposed into the contributions of the DMF chemical groups, and on the polyacrylamide chemical groups. In the first panel below we show the contributions of the DMF chemical groups to the distribution function.

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example2/script2.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

#### Output 

The decomposition reveals that specific interactions peaking at distances slightly smaller than 2$\AA$ exist between the polymer and the carbonyl group of DMF. Thus, there hydrogen bonds between the polymer and this group, which dominate the interactions between the solute and the solvent at short distances. The non-specific interactions peak at 2.5Angs and are composed of contributions of all DMF chemical groups, but particularly of the methyl groups.

![](./assets/scripts/example2/mddf_groups.png)

The decomposition of the same MDDF in the contributions of the chemical groups of the polymer is clearly associated to the DMF contributions. The specific, hydrogen-bonding, interactions, are associated to the polymer amine groups. The amine groups also contribute to the non-specific interactions at greater distances, but these are a sum of the contributions of all polymer groups, polar or aliphatic.

## [2D density map](@id 2Dmap-example2)

We can decompose the MDDF into the contributions of each portion of the polymer chain. The map below displays the contributions of each chemical group of the polymer, now split into the mers of the polymer, to the MDDF.

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example2/script3.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

#### Output 

The terminal methyl groups interact strongly with DMF, and strong local density augmentations are visible in particular on the amine groups. These occur at less than 2.0Angs and are characteristic of hydrogen-bond interactions. Interestingly, the DMF molecules are excluded from the aliphatic and carbonyl groups of the polymer, relative to the other groups.

Finally, it is noticeable that the central mer is more weakly solvated by DMF than the mers approaching the extremes of the polymer chain. This is likely a result of the partial folding of the polymer, that protects that central mers from the solvent in a fraction of the polymer configurations.

```@raw html
<center>
<img width=70% src="../assets/scripts/example2/map2D_acr.png">
</center>
```

### References

Molecules built with JSME: B. Bienfait and P. Ertl, JSME: a free molecule editor in JavaScript, Journal of Cheminformatics 5:24 (2013)
[http://biomodel.uah.es/en/DIY/JSME/draw.en.htm](http://biomodel.uah.es/en/DIY/JSME/draw.en.htm)

The system was built with [Packmol](http://m3g.iqm.unicamp.br/packmol).

The simulations were perfomed with [NAMD](https://www.ks.uiuc.edu/Research/namd/), with [CHARMM36](https://www.charmm.org) parameters. 