# [Atomic and group contributions](@id contributions)

One of the interesting features of Minimum-Distance distributions is
that they can be naturally decomposed into the atomic or group
contributions. Simply put, if a MDDF has a peak at a hydrogen-bonding
distance, it is natural to decompose that peak into the contributions of
each type of solute or solvent atom to that peak.     

To obtain the atomic contributions of an atom or group of atoms to the MDDF,
the coordination number, or the site count at each distance, the
`contributions` function is provided. For example, in a system composed
of a protein and water, we would have defined the solute and solvent
using:

```julia
using PDBTools, ComplexMixtures
atoms = readPDB("system.pdb")
protein = select(atoms,"protein")
water = select(atoms,"water")
solute = AtomSelection(protein,nmols=1)
solvent = AtomSelection(water,natomspermol=3)
```

The MDDF calculation is executed with:
```julia
trajectory = Trajectory("trajectory.dcd",solute,solvent)
results = mddf(trajectory, Options(bulk_range=(8.0, 12.0)))
```

## Atomic contributions in the result data structure

The `results` data structure contains the decomposition of the MDDF into
the contributions of every type of atom of the solute and the solvent.
These contributions can be retrieved using the `contributions` function,
with the `SoluteGroup` and `SolventGroup` selectors.

For example, if the MDDF of water (solvent) relative to a solute was computed,
and water has atom names `OH2, H1, H2`, one can retrieve the contributions
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
Each of these calls will return a vector of the constributions of these
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

- `type=:mddf` : the contribution of the group to the MDDF is returned.
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
atoms = readPDB("system.pdb")
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

## Reference functions

```@autodocs
Modules = [ComplexMixtures]
Pages = ["contributions.jl"]
```
