# [Coordination numbers](@id Tools)

- [Computing radial distribution functions](@ref radial_distribution)
- [Overview of the solvent and solute properties](@ref overview)

```@meta
CollapsedDocStrings = true
```

## [Computing radial distribution functions](@id radial_distribution)

The distributions returned by the `mddf` function (the `mddf` and
`proximal` vectors, the latter being the proximal distribution function ([Mehrotra, Marchese, Beveridge, *J. Am. Chem. Soc.* 1981](https://doi.org/10.1021/ja00393a030); [Lounnas, Pettitt, *Proteins* 1994](https://doi.org/10.1002/prot.340180206))
of the reference atom of the solvent), are normalized by the random reference state. This
is equivalent to  a site count based on the numerical integration of the volume
corresponding to each minimum-distance to the solute. 

If, however, the
solute is defined by a single atom (as the oxygen atom of water, for
example), the numerical integration of the volume can be replaced by a
simple analytical spherical shell volume, reducing noise. The `ComplexMixtures.gr`
function returns the radial distribution function and the KB integral 
computed from the results, using this volume estimate: 

```julia
g, kb = ComplexMixtures.gr(R)
```

!!! note
    The spherical shell volume is only valid if the solute has a single atom per molecule
    (the solvent may have more atoms, in which case the distances are counted relative to its
    reference atom). If the solute has more than one atom per molecule, `gr(R)` throws an 
    `ArgumentError`. The proximal distribution function ([Mehrotra, Marchese, Beveridge, *J. Am. Chem. Soc.* 1981](https://doi.org/10.1021/ja00393a030); [Lounnas, Pettitt, *Proteins* 1994](https://doi.org/10.1002/prot.340180206)) of the reference 
    atom of the solvent, normalized by the random reference state, and the corresponding KB integral, are available in the `proximal`
    and `kb_proximal` fields of the `MDDFResult` structure.

```@autodocs
Modules = [ComplexMixtures]
Pages = ["gr.jl"]
```

## [Overview of the solvent and solute properties](@id overview)

The output to the REPL of the MDDFResult structure provides an overview of the
properties of the solution. The data can be retrieved into a data structure
using the `overview` function. Examples:     

```julia-repl
...

julia> results = mddf(trajectory_file, solute, solvent, Options(cutoff=12.0))

julia> results
--------------------------------------------------------------------------------
MDDF Overview - ComplexMixtures - Version 2.0.8
--------------------------------------------------------------------------------

Solvent properties:
-------------------

Simulation concentration: 0.49837225882780106 mol L⁻¹
Molar volume: 2006.532230249041 cm³ mol⁻¹

Concentration in bulk: 0.5182380507741433 mol L⁻¹
Molar volume in bulk: 1929.6151614228274 cm³ mol⁻¹

Solute properties:
------------------

Simulation Concentration: 0.002753437894076249 mol L⁻¹
Estimated solute partial molar volume: 13921.98945754469 cm³ mol⁻¹

Bulk range: 8.0 - 12.0 Å
Molar volume of the solute domain: 34753.1382279134 cm³ mol⁻¹

Auto-correlation: false

Trajectory files and weights:

   /home/user/NAMD/trajectory.dcd - w = 1.0

Long range MDDF mean (expected 1.0): 1.0378896753018338 ± 1.0920172247127446
Long range RDF mean (expected 1.0): 1.2147429551790854 ± 1.2081838161780682

--------------------------------------------------------------------------------
```

In this case, since solute and solvent are equivalent and the system is
homogeneous, the molar volumes and concentrations are similar. This is
not the case if the molecules are different or if the solute is at
infinite dilution (in which case the bulk solvent density might be
different from the solvent density in the simulation). 

To retrieve the data of the overview structure use, for example:

```julia-repl
julia> overview = overview(results);

julia> overview.solute_molar_volume
657.5051512801567
```

## [Renormalization of bulk density](@id renormalize)

The `renormalize` function allows to renormalize the `MDDFResult` structure to a different bulk density of the solvent. This is useful in cases where the bulk density cannot be estimated from the simulation, for example in simulations of zeolites, MOFs, etc, where the solvent is confined within a porous material and there's no bulk region.

Under the hood, the function rescales the site count of the reference (ideal gas) distribution used for normalization.

Typically, the MDDF would be computed with an arbitrary bulk region and renormalized using the density of the solvent in the complete simulation box, or the density of pure solvent.

```@docs
renormalize
```

!!! warning 
    In these situations the computation of KB integrals, apparent molar volumes, and preferential interactions might not be meaningful, but the MDDF and RDF can still provide useful information about the local structure of the solvent around the solute.




