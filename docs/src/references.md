# References

## Primary citations

If this package was useful to you, please cite the following papers:

* L. Martínez, **ComplexMixtures.jl: Investigating the structure of solutions of complex-shaped molecules from a solvent-shell perspective.** *J. Mol. Liq.* 347, 117945, 2022.  [[Full Text]](https://doi.org/10.1016/j.molliq.2021.117945)

* L. Martínez, S. Shimizu, **Molecular interpretation of preferential interactions in protein solvation: a solvent-shell perspective by means of minimum-distance distribution functions.** *J. Chem. Theor. Comp.* 13, 6358–6372, 2017. [[Full Text]](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00599)

## Applications and examples

* A. F. Pereira, L. Martínez, **Osmolyte Structural and Thermodynamic Effects Across the Protein Folding Landscape**, *JACS Au*, 2025. [[Full Text]](https://doi.org/10.1021/jacsau.5c00813). 

```@raw html
<a href="https://pubs.acs.org/toc/jpcbfk/129/27"><img style="float: right" src="../assets/2025_jpcb_capa.png" width="20%"></a>
```

* V. Piccoli, L. Martínez, **Cation Hydrophobicity Effects on Protein Solvation in Aqueous Ionic Liquids**, *J. Phys. Chem. B*, 129 (27) 6765-6776, 2025. [[Full Text]](https://pubs.acs.org/doi/10.1021/acs.jpcb.5c00779). 

* F. C. Ramos, L. Martínez, **Molecular dynamics and solvation structures of the β-glucosidase from Humicola insolens (BGHI) in aqueous solutions containing glucose**
  *Int. J. Biol. Macromol.* 286 (138210) 2025. [[Full Text]](https://doi.org/10.1016/j.ijbiomac.2024.138210)

* V. Piccoli, L. Martínez, **Competitive Effects of Anions on Protein Solvation by Aqueous Ionic Liquids.** 
  *J. Phys. Chem. B* 128, 7792-7802, 2024. [[Full Text]](https://pubs.acs.org/doi/10.1021/acs.jpcb.4c03735)

* A. F. Pereira, L. Martínez, **Helical Content Correlations and Hydration Structures of the Folding Ensemble of the B Domain of Protein A.**
  *J. Chem. Inf. Model.* 64, 3350-3359, 2024. [[Full Text]](https://doi.org/10.1021/acs.jcim.3c01822)

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

## See also

* [Packmol](http://m3g.iqm.unicamp.br/packmol): A package for building initial configurations for molecular dynamics simulations.

* [CellListMap.jl](https://m3g.github.io/CellListMap.jl/stable/): Efficient and customizable implementation of cell lists, which allows the computation of general properties dependent on distances of particles within a cutoff, for example short-range potentials, forces, neighbor lists, etc.

* [MDLovoFit](http://m3g.iqm.unicamp.br/mdlovofit/home.shtml): Automatic identification of mobile and rigid substructures in molecular dynamics simulations and fractional structural fluctuation analysis. 

## Breaking changes

The syntax changes necessary to update script from version `1.X` to `2.X` of 
the package are:

### Atom selections

The previous `Selection` structure was renamed to `AtomSelection` for clarity.
- Before:
```julia
water = Selection(water; natomspermol=3)
```
- Now:
```julia
water = AtomSelection(water; natomspermol=3)
```

### Group contributions syntax

The syntax to computing group contributions is improved. Previously, the `contrib` or
`contributions` functions required three somewhat redundant parameters. 
- Before:
The call to `contributions` required 3 parameters: the `Selection` structure,
the matrix of contributions, and the indexes of the atoms for which the
contributions were desired:
```julia
h_contributions = contributions(solvent, R.solvent_atom, h_indexes)
```
- Now:
The contributions are extracted from the `Result` data structure, by 
providing either a `SoluteGroup` or `SolventGroup` object, which are
setup with the group names, group indexes, atom names, or atom indexes:
```julia
h_contributions = contributions(R, SolventGroup(h_indexes))
```

### Frame weights

`frame_weights` is an option of the `mddf` execution. That is previously,
they were defined in the `Options` data structure, and now they are passed
to the `mddf` function.
- Before:
```julia
options = Options(frame_weights=[1.0, 2.0], bulk_range=(8.0, 12.0))
results = mddf(trajectory_file, solute, solvent, options)
```
- Now:
```julia
results = mddf(trajectory_file, solute, solvent, options; frame_weights=[1.0, 2.0])
```
