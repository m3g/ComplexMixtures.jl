# Introduction

[ComplexMixtures.jl](https://github.com/m3g/ComplexMixtures.jl) is a package to study the solute and solvent interactions of
mixtures of molecules of complex shape. Conventional radial distribution
functions are not appropriate to represent the structure of a solvent
around a solute with many atoms, and a variable, non-spherical shape.     

Typical solutes of complex shape are proteins, nucleic acids, and
polymers in general. Smaller molecules like lipids, carbohydrates, etc,
are also complex enough such that representing the structure of the
solution of those molecules with distribution functions is not trivial.

Minimum-Distance Distribution Functions (MDDFs) are a very general and
practical way to represent solute-solvent interactions for molecules
with arbitrarily complex sizes and geometries. Briefly, instead of
computing the density distribution function of a particular atom or the
center-of-mass of the molecules, one computes the distribution function
of the minimum-distance between any solute and solvent atoms. This
provides a size and shape-independent distribution which is very natural
to interpret in terms of molecular interactions.   

Additionally, the MDDFs can be decomposed into contributions of each
type of atom (or groups of atoms) of the solute and solvent molecules,
such that the profiles of the distributions can be interpreted in terms
of the chemical nature of the species involved in the interactions at
each distance.   

Finally, as with radial distribution functions, MDDFs can be used to
compute Kirkwood-Buff integrals to connect the accumulation or depletion
of the solvents components to thermodynamic properties, like protein
structural stability, solubility, and others.

## Features

Check out our [examples repository](https://github.com/m3g/ComplexMixturesExamples), featuring the analysis of solvation structures for proteins, polymers, membrane, and complex solutions! The examples are also described in our [featured article](https://doi.org/10.1016/j.molliq.2021.117945).

### 1. Minimum-distance distribution functions: understanding solvation at a molecular level

This figure illustrates one of the main features of minimum-distance distribution functions, by showing the distribution of DMF molecules at the surface of an polyacrylamide molecule. The direct interactions are evident by the peak at hydrogen-bonding distances and, additionally, the contribution of each group of atoms of the DMF can be clearly distinguished by decomposing the total MDDF into atomic or chemical group contributions. 

```@raw html
<center>
<img width=60% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Polyacrylamide_in_DMF/results/mddf_groups.png">
<br>
<b>Minimum distance distribution function and its decomposition into the chemical
groups of the solvent (top) and solute (bottom) molecules.<br><br></b> 
</center>
```

Decomposition of the total MDDF into the contributions of the solute atoms (in this case, a protein) is also possible. Any chemical group decomposition is possible. Here, we decompose the MDDF into the contribution of each protein residue. 

```@raw html
<center>
<img width=70% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/Density2D/density2D.png">
<br>
<b>Density map of a solvent in the vicinity of each protein residue.</b> 
</center>
```

### 2. Thermodynamic interpretation through Kirkwood-Buff theory

Minimum-distance distribution functions can be used to compute Kirkwood-Buff integrals, and thus, thermodynamic parameters associated to solvation. 

Kirkwood-Buff integrals carry the information of the total accumulation or depletion of each solvent around a solute. For example, the figure below displays the KB integrals of an ionic liquid solvating different conformational states of a protein [[link]](https://www.sciencedirect.com/science/article/pii/S016773222201491X?via%3Dihub). The figure illustrates that the solvation structures are dependent on the protein folding state. 

```@raw html
<center>
<img width="50%" src="./figures/KBs.png">
<br>
<b>Kirkwood-Buff integrals of an ionic liquid solvating a protein in different conformational states.</b><br><br> 
</center>
```

From differences in KB integrals among cossolvents, the Preferential Solvation parameter can be computed. This is an important parameter because it can be measured experimentally and is ultimately associated with the equilibrium thermodynamics of the solvation. In the following figure, we show that, for example, the preferential solvation of a protein in different folding states is dependent in a non-trivial way on the concentration of an ionic liquid in aqueous solutions. 

```@raw html
<center>
<img width="50%" src="./figures/Gamma.png">
<br>
<b>Kirkwood-Buff integrals of an ionic liquid solvating a protein in different conformational states.</b><br><br>
</center>
```

In particular, the plot shows that besides being preferentially excluded from the protein surface at high concentrations in the native state, suggesting protein folding stabilization, the interactions with the protein in the denatured states are stronger, leading to denaturation at all concentrations. 


## References

* L. Martínez, **ComplexMixtures.jl: Investigating the structure of solutions of complex-shaped molecules from a solvent-shell perspective.** *J. Mol. Liq.* 117945, 2021.  [[Full Text]](https://doi.org/10.1016/j.molliq.2021.117945)

* L. Martínez, S. Shimizu, **Molecular interpretation of preferential interactions in protein solvation: a solvent-shell perspective by means of minimum-distance distribution functions.** *J. Chem. Theor.  Comp.* 13, 6358–6372, 2017. [[Full Text]](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00599)

## See also

### Seminar

* Presentation about ComplexMixtures.jl and protein-solvent interactions: [https://youtu.be/umSRjsITzyA](https://youtu.be/umSRjsITzyA)

### Applications

* A. F. Pereira, V. Piccoli, L. Martínez, **Trifluoroethanol direct interactions with protein backbones destabilize alpha-helices.** 
  *J. Mol. Liq.* 365, 120209, 2022. [[Full Text]](https://doi.org/10.1016/j.molliq.2022.120209)

* V. Piccoli, L. Martínez, **Ionic liquid solvation of proteins in native and denatured states.** 
  *J. Mol. Liq.* 363, 119953, 2022. [[Full Text]](http://dx.doi.org/10.1016/j.molliq.2022.119953)

* V. Piccoli, L. Martínez, **Correlated counterion effects in the solvation of proteins by ionic-liquids.** *J. Mol. Liq.* 320, 114347, 2020.
  [[Full Text]](https://doi.org/10.1016/j.molliq.2020.114347)

* I. P. de Oliveira, L. Martínez, **The shift in urea orientation at protein surfaces at low pH is compatible with a direct mechanism of protein denaturation.** *Phys. Chem. Chem. Phys.* 22, 354-367, 2020.
  [[Full Text]](https://pubs.rsc.org/en/content/articlelanding/2019/CP/C9CP05196A#!divAbstract)

* I. P. de Oliveira, L. Martínez, **Molecular basis for competitive solvation of the Burkholderia cepacia lipase by sorbitol and urea.**
  *Phys. Chem. Chem. Phys.* 18, 21797-21808, 2016.
  [[Full Text]](https://pubs.rsc.org/en/content/articlelanding/2016/cp/c6cp01789d#!divAbstract)

