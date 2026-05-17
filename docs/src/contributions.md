```@meta
CollapsedDocStrings = true
```

# [Atomic and group contributions](@id contributions)

One of the interesting features of Minimum-Distance distributions is
that they can be naturally decomposed into the atomic or group
contributions. Simply put, if a MDDF has a peak at a hydrogen-bonding
distance, it is natural to decompose that peak into the contributions of
each type of solute or solvent atom to that peak.     

!!! tip
    See also the section on contributions per residue for proteins and other
    macromolecules: [2D density map per residue](@ref 2D_per_residue).

To obtain the atomic contributions of an atom or group of atoms to the MDDF,
the coordination number, or the site count at each distance, the
`contributions` function is provided. For example, in a system composed
of a protein and water, we would have defined the solute and solvent
using:

```julia
using PDBTools, ComplexMixtures
atoms = read_pdb("system.pdb")
protein = select(atoms,"protein")
water = select(atoms,"water")
solute = AtomSelection(protein,nmols=1)
solvent = AtomSelection(water,natomspermol=3)
```

The MDDF calculation is executed with:
```julia
results = mddf("trajectory.dcd", solute, solvent, Options(bulk_range=(8.0, 12.0)))
```

## Atomic contributions in the result data structure

The `results` data structure contains the decomposition of the MDDF into
the contributions of every type of atom of the solute and the solvent.
These contributions can be retrieved using the `contributions` function,
with the `SoluteGroup` and `SolventGroup` selectors.

```@docs
contributions
SolventGroup
SoluteGroup
```

### Example: computing the oxygen contributions of water

Here we show the MDDF of water (solvent) relative to a solute.
Water molecules have atom names `OH2, H1, H2`, one can retrieve the contributions
of the oxygen atom with:

```julia
OH2 = contributions(results, SolventGroup(["OH2"]))
```
or with, if `OH2` is the first atom in the molecule,
```julia
OH2 = contributions(results, SolventGroup([1]))
```

The contributions of the hydrogen atoms can be obtained, similarly, with:
```julia
H = contributions(results, SolventGroup(["H1", "H2"]))
```
or with, if `OH2` is the first atom in the molecule,
```julia
H = contributions(results, SolventGroup([2, 3]))
```
Each of these calls will return a vector of the contributions of these
atoms to the total MDDF. 

For example, here we plot the total MDDF and the Oxygen contributions: 

```julia
using Plots
plot(results.d, results.mddf, label=["Total MDDF"], linewidth=2)
plot!(results.d, contributions(results, SolventGroup(["OH2"])), label=["OH2"], linewidth=2)
plot!(xlabel="Distance / Å", ylabel="MDDF")
```

```@raw html
<img src="../figures/oh2.png" width="60%">
```
## Contributions to coordination numbers or site counts

The keyword `type` defines the return type of the contribution:

- `type=:mddf` : the contribution of the group to the MDDF is returned (default).
- `type=:coordination_number` : the contribution of the group to the coordination number, that is, the 
   cumulative sum of counts at each distance, is returned.
- `type=:md_count` : the contribution of the group to the site count at each distance is returned. 

Example of the usage of the `type` option:
```julia
ca_contributions = contributions(results, SoluteGroup(["CA"]); type=:coordination_number)
```

## Using PDBTools

If the solute is a protein, or other complex molecule, selections defined
with `PDBTools` can be used. For example, this will retrieve the contribution
of the acidic residues of a protein to total MDDF:
```julia
using PDBTools
atoms = read_pdb("system.pdb")
acidic_residues = select(atoms, "acidic")
acidic_contributions = contributions(results, SoluteGroup(acidic_residues))
```
It is expected that for a protein most of the atoms do not contribute to
the MDDF, and that all values are zero at very short distances, smaller
than the radii of the atoms.

More interesting and general is to select atoms of a complex
molecule, like a protein, using residue names, types, etc. Here we
illustrate how this is done by providing selection strings to
`contributions` to obtain the contributions to the MDDF of different
types of residues of a protein to the total MDDF. 

For example, if we want to split the contributions of the charged and
neutral residues to the total MDDF distribution, we could use to following
code. Here, `solute` refers to the protein.

```julia
charged_residues = PDBTools.select(atoms,"charged")
charged_contributions = contributions(results, SoluteGroup(charged_residues))

neutral_residues = PDBTools.select(atoms,"neutral")
neutral_contributions = contributions(atoms, SoluteGroup(neutral_residues))
```

The `charged_contributions` and `neutral_contributions` outputs are vectors containing the
contributions of these residues to the total MDDF. The corresponding
plot is:   

```julia
plot(results.d,results.mddf,label="Total MDDF",linewidth=2)
plot!(results.d,charged_contributions,label="Charged residues",linewidth=2)
plot!(results.d,neutral_contributions,label="Neutral residues",linewidth=2)
plot!(xlabel="Distance / Å",ylabel="MDDF")
```
Resulting in:

```@raw html
<img src="../figures/charged_and_neutral.png" width="60%">
```

Note here how charged residues contribute strongly to the peak at
hydrogen-bonding distances, but much less in general. Of course all
selection options could be used, to obtain the contributions of specific
types of residues, atoms, the backbone, the side-chains, etc. 

## Proximal contributions to KBIs

The `type=:kbi` option in `contributions` computes the **proximal contribution** of each atomic group to the
Kirkwood-Buff integral (KBI). Each solute atom acts as a reference: solvent molecules for which
that atom is the nearest solute atom are assigned to it. The KBI contribution of a group is then computed
from the excess or deficit of proximal solvent contacts relative to the random (ideal-gas) reference,
normalized by the bulk solvent density — exactly as the total KBI is computed, but restricted to the
subset of contacts attributed to that group. The result is a decomposition of the total KBI (`results.kb`)
into atomic or group contributions:

```math
G_{total}(r) = \sum_i G_i(r)
```

where each $G_i(r)$ is the proximal contribution of group $i$.

For example, consider a protein solvated by Glucose (`resname GLYC`), and follow the  `mddf` computation steps of [this example](@ref example1): 

```julia
using ComplexMixtures, PDBTools, Plots, LaTeXStrings
atoms = read_pdb("system.pdb")
solute = AtomSelection(select(atoms, "protein"); nmols=1)
solvent = AtomSelection(select(atoms, "resname GLYC"); natomspermol=14)
R = mddf("glyc50_traj.dcd", solute, solvent, Options(bulk_range=(10,12)))
```

The contributions of charged, polar but not charged, and non-polar residues to the KBI can be computed with:

```julia
charged_atoms = select(atoms, "protein and charged")
polar_not_charged_atoms = select(atoms, "protein and polar and not charged")
nonpolar_atoms = select(atoms, "protein and nonpolar")
kbi_charged = contributions(R, SoluteGroup(charged_atoms); type=:kbi)
kbi_polar_not_charged = contributions(R, SoluteGroup(polar_not_charged_atoms); type=:kbi)
kbi_nonpolar = contributions(R, SoluteGroup(nonpolar_atoms); type=:kbi)
```

These partial contributions are additive: `R.kb ≈ kbi_charged + kbi_polar_not_charged + kbi_nonpolar`.

```julia
plot(R.d, R.kb / 1000, label="Total KBI", linewidth=2)
plot!(R.d, kbi_charged / 1000, label="Charged residues", linewidth=2)
plot!(R.d, kbi_polar_not_charged / 1000, label="Polar, not charged, residues", linewidth=2)
plot!(R.d, kbi_nonpolar / 1000, label="Non-polar residues", linewidth=2)
plot!(xlabel="Distance / Å", ylabel="KBI / L mol⁻¹")
```

```@raw html
<img src="../figures/kbi_contributions1.png" width="60%">
```

!!! note
    The KBI contributions are returned in cm³ mol⁻¹, consistent with `results.kb`. Divide by 1000 to convert to L mol⁻¹.

## Per-residue proximal contributions to KBIs

The `ResidueContributions` function supports `type=:kbi`, producing a 2D map of the proximal KBI
contributions decomposed per residue. Since the KBI requires integration over distances long enough to
converge to the bulk, the `dmax` parameter should be extended well beyond its default value of 3.5 Å. Continuing from the previous example, we have:

```julia
rc_kbi = ResidueContributions(R, protein; type=:kbi, dmax=12.0)
heatmap(rc_kbi)
```

```@raw html
<img src="../figures/kbi_contributions2.png" width="80%">
```

The sum of contributions from all residues at the last distance converges to the total KBI:

```julia
sum(rc_kbi[i][end] for i in eachindex(rc_kbi)) ≈ results.kb[end]
```

Or, since the converged value of the KBIs is of particular interest, a bar plot of the proximal contributions at the final distance can illustrate better the contribution of each residue:

```julia
bar(
    [rc_kbi[i][end]/1000 for i in eachindex(rc_kbi)]; 
    xticks=PDBTools.residue_ticks(protein; stride=10), xrotation=60, 
    label="", 
    framestyle=:box, 
    ylabel=L"G_i~/~\textrm{L~mol^{-1}}", 
    xlabel="residue"
)
```

```@raw html
<img src="../figures/kbi_contributions3.png" width="80%">
```

Residues with positive contributions are those for which local accumulation of the solvent compensates the excluded volume effect. Residues that are completely hidden from the surface will contribute proportionally to their excluded volumes.

!!! tip
    See the [2D density map per residue](@ref 2D_per_residue) section for details on indexing, slicing, and arithmetic operations on `ResidueContributions` objects.

