# ComplexMixtures

A package to study the structure of solutions formed by solutes and solvents of complex molecular shapes.

## Documentation:

The documentation is available at: [http://m3g.iqm.unicamp.br/ComplexMixtures](http://m3g.iqm.unicamp.br/ComplexMixtures)

## Overview

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

<center>
<img width=70% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/MDDF/mddf_atom_contrib.png">
<br>
<b>Minimum-distance distribution function and its decomposition into molecular groups.</b> 
</center>

Additionally, the MDDFs can be decomposed into contributions of each
type of atom (or groups of atoms) of the solute and solvent molecules,
such that the profiles of the distributions can be interpreted in terms
of the chemical nature of the species involved in the interactions at
each distance.   

Finally, as with radial distribution functions, MDDFs can be used to
compute Kirkwood-Buff integrals to connect the accumulation or depletion
of the solvents components to thermodynamic properties, like protein
structural stability, solubility, and others.

<center>
<img width=70% src="https://github.com/m3g/ComplexMixturesExamples/raw/main/Density2D/density2D.png">
<br>
<b>Density map of a solvent in the vicinity of each protein residue..</b> 
</center>

### References

* L. Martínez, **ComplexMixtures.jl: Understanding the solute-solvent interactions
  of molecules of complex molecular shapes.** To be published, 2021.

* L. Martínez, S. Shimizu, **Molecular interpretation of preferential
  interactions in protein solvation: a solvent-shell perspective by
  means of minimum-distance distribution functions.** *J. Chem. Theor.
  Comp.* 13, 6358–6372, 2017. [[Full Text]](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00599)




