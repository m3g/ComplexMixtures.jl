# Introduction

*ComplexMixtures* is a package to study the solute and solvent interactions of
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

```@raw html
<center>
<img width=70% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Protein_in_Glycerol/Density2D/density2D.png">
<br>
<b>Density map of a solvent in the vicinity of each protein residue (see the Full Example section).</b> 
</center>
```

## References

* L. Martínez, **ComplexMixtures.jl: Investigating the structure of solutions of complex-shaped molecules from a solvent-shell perspective.** *J. Mol. Liq.* 117945, 2021.  [[Full Text]](https://doi.org/10.1016/j.molliq.2021.117945)

* L. Martínez, S. Shimizu, **Molecular interpretation of preferential interactions in protein solvation: a solvent-shell perspective by means of minimum-distance distribution functions.** *J. Chem. Theor.  Comp.* 13, 6358–6372, 2017. [[Full Text]](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00599)

## See also

### Seminar

* Presentation about ComplexMixtures.jl and protein-solvent interactions: [https://youtu.be/umSRjsITzyA](https://youtu.be/umSRjsITzyA)

### Applications

* A. F. Pereira, V. Piccoli, L. Martínez, **Trifluoroethanol direct interactions with protein backbones destabilize alpha-helices.** *J. Mol. Liq.* 365 (2022) 120209. [[Full Text]](https://doi.org/10.1016/j.molliq.2022.120209)

* V. Piccoli, L. Martínez, **Ionic liquid solvation of proteins in native and denatured states.** *J. Mol. Liq.* 363 (2022) 119953. [[Full Text]](http://dx.doi.org/10.1016/j.molliq.2022.119953)

* V. Piccoli, L. Martínez, **Correlated counterion effects on the solvation of proteins by ionic liquids.** *J. Mol. Liq.* 320, 114347, 2020. [[Full Text]](https://www.sciencedirect.com/science/article/pii/S0167732220337247?dgcid=author)

* I. P. de Oliveira, L. Martínez, **The shift in urea orientation at protein surfaces at low pH is compatible with a direct mechanism of protein denaturation.** *Phys. Chem. Chem. Phys.* 22, 354-367, 2020. [[Full Text]](https://pubs.rsc.org/en/content/articlelanding/2019/CP/C9CP05196A#!divAbstract)
 
