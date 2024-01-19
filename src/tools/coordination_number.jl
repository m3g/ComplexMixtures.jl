"""
    coordination_number(R::Result) = R.coordination_number
    coordination_number(R::Result, group_contributions::Vector{Float64})
    coordination_number(s::AtomSelection, atom_contributions::Matrix{Float64}, R::Result, group)

Computes the coordination number of a given group of atoms from the solute or solvent atomic contributions to the MDDF. 
If no group is defined (first call above), the coordination number of the whole solute or solvent is returned.

If the `group_contributions` to the `mddf` are computed previously with the `contributions` function, the result can be used
to compute the coordination number by calling `coordination_number(R::Result, group_contributions)`.

Otherwise, the coordination number can be computed directly with the second call, where:

`s` is the solute or solvent selection (type `ComplexMixtures.AtomSelection`)

`atom_contributions` is the `R.solute_atom` or `R.solvent_atom` arrays of the `Result` structure

`R` is the `Result` structure,

and the last argument is the selection of atoms from the solute to be considered, given as a list of indices, list of atom names, 
or a selection following the syntax of `PDBTools`, or vector of `PDBTools.Atom`s, or a `PDBTools.Residue`

# Examples

In the following example we compute the coordination number of the atoms of residue 50 (of the solute) with the solvent atoms of TMAO,
as a function of the distance. Finally, we show the average number of TMAO molecules within 5 Angstroms of residue 50. 
The `findlast(<(5), R.d)` part of the code below returns the index of the last element of the `R.d` array that is smaller than 5 Angstroms.

## Precomputing the group contributions Using the `contributions` function

```julia
using ComplexMixtures, PDBTools
pdb = readPDB("test/data/NAMD/structure.pdb");
R = load("test/data/NAMD/protein_tmao.json");
solute = AtomSelection(PDBTools.select(pdb, "protein"), nmols=1);
residue50 = PDBTools.select(pdb, "residue 50");
# Compute the group contributions to the MDDF
residue50_contribution = contributions(solute, R.solute_atom, residue50);
# Now compute the coordination number
residue50_coordination = coordination_number(R, residue50_contribution)
# Output the average number of TMAO molecules within 5 Angstroms of residue 50
residue50_coordination[findlast(<(5), R.d)]
```

## Without precomputing the `group_contribution`

```julia
using ComplexMixtures, PDBTools
pdb = readPDB("test/data/NAMD/structure.pdb");
R = load("test/data/NAMD/protein_tmao.json");
solute = AtomSelection(PDBTools.select(pdb, "protein"), nmols=1);
residue50 = PDBTools.select(pdb, "residue 50");
# Compute the coordination number
residue50_coordination = coordination_number(solute, R.solute_atom, R, group)
# Output the average number of TMAO molecules within 5 Angstroms of residue 50
residue50_coordination[findlast(<(5), R.d)]
```

"""
function coordination_number end

coordination_number(R::Result) = R.coordination_number

function coordination_number(
    s::AtomSelection,
    atom_contributions::Matrix{Float64},
    R::Result,
    group,
)
    # Extract the group contributions to the MDDF
    group_contributions = contributions(s, atom_contributions, group)
    # Compute the coordination number
    return coordination_number(R, group_contributions)
end
function coordination_number(R::Result, group_contributions::Vector{Float64})
    # obs: we accumulate group coordination numbers, but at the end we 
    # divide by random count of, so we the group contributions are the contributions 
    # to the MDDF. To recover the coordination number we need to multiply back
    # the group contributions by the random count.
    cn = cumsum(
        group_contributions[i] * R.md_count_random[i] for
        i in eachindex(group_contributions)
    )
    return cn
end

@testitem "coordination_number" begin
    using ComplexMixtures, PDBTools
    using ComplexMixtures.Testing
    dir = "$(Testing.data_dir)/NAMD"
    atoms = readPDB("$dir/structure.pdb")
    protein = AtomSelection(select(atoms, "protein"), nmols = 1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol = 14)
    options = Options(lastframe = 1, seed = 321, StableRNG = true, nthreads = 1, silent = true, n_random_samples=200)
    traj = Trajectory("$dir/trajectory.dcd", protein, tmao)
    R = mddf(traj, options)

    # Test self-consistency
    @test sum(R.solute_atom) ≈ sum(R.solvent_atom)
    @test coordination_number(R, contributions(protein, R.solute_atom, select(atoms, "protein"))) == R.coordination_number

    # Checked with vmd: same residue as (resname TMAO and within 3.0 of protein)
    @test R.coordination_number[findfirst(>(3), R.d)] == 7.0
    @test R.coordination_number[findfirst(>(5), R.d)] == 14.0

    # THR93 forms a hydrogen bond with TMAO in this frame
    thr93 = select(atoms, "protein and residue 93")
    @test last(coordination_number(R, contributions(protein, R.solute_atom, thr93))) == 1

    # Consistency of a read-out result
    pdb = readPDB("$(Testing.data_dir)/NAMD/Protein_in_Glycerol/system.pdb")
    R = load("$(Testing.data_dir)/NAMD/Protein_in_Glycerol/protein_water.json"; legacy_warning=false)
    group = PDBTools.select(pdb, "protein")
    solute = AtomSelection(group, nmols = 1)
    cn = coordination_number(solute, R.solute_atom, R, group)
    @test maximum(R.coordination_number .- cn) < 1e-10
end
