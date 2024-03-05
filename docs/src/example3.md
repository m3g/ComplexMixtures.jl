# POPC membrane in water/ethanol

In this example ComplexMixtures.jl is used to study the interactions of a POPC membrane with a mixture of 20%(mol/mol) ethanol in water. At this concentration ethanol destabilizes the membrane.

```@raw html
<center>
<img width=50% src="https://raw.githubusercontent.com/m3g/ComplexMixturesExamples/main/POPC_in_Water-Ethanol/system.png">
</center>
```

System image: a POPC membrane (center) solvated by a mixture of water (purple) and ethanol (green). The system is composed by 59 POPC, 5000 water, and 1000 ethanol molecules.  

### Index

- [Data, packages, and execution](@ref data-example3)
- [MDDF and KB integrals](@ref mddf-example3)
- [Group contributions](@ref groups1-example3)
- [Interaction of POPC groups with water](@ref groups2-example3)
- [Interaction of POPC groups with ethanol](@ref groups3-example3)
- [Density map on POPC chains](@ref map-example3)

## [Data, packages, and execution](@id data-example3)

The files required to run this example are:

- [equilibrated.pdb](https://raw.githubusercontent.com/m3g/ComplexMixturesExamples/main/POPC_in_Water-Ethanol/simulation/equilibrated.pdb): The PDB file of the complete system.
- [traj_POPC.dcd](https://www.dropbox.com/scl/fi/hcenxrdf8g8hfbllyakhy/traj_POPC.dcd?rlkey=h9zivtwgya3ivva1i6q6xmr2p&dl=0): Trajectory file. This is a 365Mb file, necessary for running from scratch the calculations.

To run the scripts, we suggest the following procedure:

1. Create a directory, for example `example3`.
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

## [MDDF and KB integrals](@id mddf-example3)

Here we show the distribution functions and KB integrals associated to the solvation of the membrane by water and ethanol. 

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example3/script1.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

#### Output 

The distribution functions are shown in the first panel of the figure below, and the KB integrals are shown in the second panel.

![](./assets/scripts/example3/popc_water_ethanol_mddf_kb.png)

Clearly, both water and ethanol accumulate on the proximity of the membrane. The distribution functions suggest that ethanol displays a greater local density augmentation, reaching concentrations roughly 4 times higher than bulk concentrations. Water has a peak at hydrogen-bonding distances (~1.8$$\mathrm{\AA}$$) and a secondary peak at 2.5$$\mathrm{\AA}$$.

Despite the fact that ethanol displays a greater relative density (relative to its own bulk concentration) at short distances, the KB integral of water turns out to be greater (more positive) than that of ethanol. This implies that the membrane is preferentially hydrated.

## [Ethanol group contributions](@id groups1-example3)

The minimum-distance distribution function can be decomposed into the contributions of the ethanol molecule groups. 

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example3/script2.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

In the figure below we show the contributions of the ethanol hydroxyl and aliphatic chain groups to the total MDDF.

![](./assets/scripts/example3/mddf_ethanol_groups.png)

As expected, the MDDF at hydrogen-bonding distances is composed by contributions of the ethanol hydroxyl group, and the non-specific interactions at ~2.5$$\mathrm{\AA}$$ have a greater contribution of the aliphatic chain of the solvent molecules. It is interesting to explore the chemical complexity of POPC in what concerns these interactions.

## [Interaction of POPC groups with water](@id groups2-example3)

The MDDF can also be decomposed into the contributions of the solute atoms and chemical groups. First, we show the contributions of the POPC chemical groups to the water-POPC distribution.

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example3/script3.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

![](./assets/scripts/example3/mddf_POPC_water_groups.png)

Not surprisingly, water interactions occur majoritarily with the Phosphate and Choline groups of POPC molecules, that is, with the polar head of the lipid. The interactions at hydrogen-bonding distances are dominated by the phosphate group, and non-specific interaction occur mostly with the choline group. Some water molecules penetrate the membrane and interact with the glycerol and aliphatic chains of POPC, but these contributions are clearly secondary.

## [Interaction of POPC groups with ethanol](@id groups3-example3)

The interactions of ethanol molecules with the membrane are more interesting, because ethanol penetrates the membrane. Here we decompose the ethanol-POPC distribution function into the contributions of the POPC chemical groups.

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example3/script4.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

![](./assets/scripts/example3/mddf_POPC_ethanol_groups.png)

Ethanol molecules interact with the choline and phosphate groups of POPC molecules, as do water molecules. The contributions to the MDDF at hydrogen-bonding distances come essentially from ethanol-phosphate interactions.

However, ethanol molecules interact frequently with the glycerol and aliphatic chains of POPC. Interactions with the Oleoyl chain are slightly stronger than with the Palmitoyl chain. This means that ethanol penetrates the hydrophobic core of the membrane, displaying non-specific interactions with the lipids and with the glycerol group. These interactions are probably associated to the destabilizing role of ethanol in the membrane structure.

## [Density map on POPC chains](@id map-example3)

The MDDFs can be decomposed at more granular level, in which each chemical group of the aliphatic chains of the POPC molecules are considered independently. This allows the study of the penetration of the ethanol molecules in the membrane. In the figure below, the carbonyl following the glycerol group of the POPC molecules is represented in the left, and going to the right the aliphatic chain groups are sequentially shown.

```@raw html
<details><summary><font color="darkgreen">Complete example code: click here!</font></summary>
```
```@eval
using Markdown
code = Markdown.parse("""
\`\`\`julia
$(read("./assets/scripts/example3/script5.jl", String))
\`\`\`
""")
```
```@raw html
</details><br>
```

![](./assets/scripts/example3/POPC_ethanol_chains.png)

Ethanol displays an important density augmentation at the vicinity of the carbonyl that follows the glycerol group, and accumulates on the proximity of the aliphatic chain. The density of ethanol decreases as one advances into the aliphatic chain, displaying a minimum around the insaturation in the Oleoyl chain. The terminal methyl group of both chains display a greater solvation by ethanol, suggesting the twisting of the aliphatic chain expose these terminal groups to membrane depth where ethanol is already abundant.

The equivalent maps for water are strikingly different, and show that water is excluded from the interior of the membrane:

![](./assets/scripts/example3/POPC_water_chains.png)

### References

Membrane built with the [VMD](https://www.ks.uiuc.edu/Research/vmd/) membrane plugin. 

Water and ethanol layers added with [Packmol](http://m3g.iqm.unicamp.br/packmol).

The simulations were performed with [NAMD](https://www.ks.uiuc.edu/Research/namd/), with [CHARMM36](https://www.charmm.org) parameters. 

Density of the ethanol-water mixture from: [https://wissen.science-and-fun.de/chemistry/chemistry/density-tables/ethanol-water-mixtures/](https://wissen.science-and-fun.de/chemistry/chemistry/density-tables/ethanol-water-mixtures/)



