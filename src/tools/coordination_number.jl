"""
    coordination_number(R::Result, group_contributions::Vector{Float64})
    coordination_number(s::Selection, atom_contributions::Matrix{Float64}, R::Result, group)

Computes the coordination number of a given group of atoms from the solute or solvent atomic contributions to the MDDF.

If the `group_contributions` to the `mddf` are computed previously with the `contrib` function, the result can be used
to compute the coordination number by calling `coordination_number(R::Result, group_contributions)`.

Otherwise, the coordination number can be computed directly with the second call, where:

`s` is the solute or solvent selection (type `ComplexMixtures.Selection`)

`atom_contributions` is the `R.solute_atom` or `R.solvent_atom` arrays of the `Result` structure

`R` is the `Result` structure,

and the last argument is the selection of atoms from the solute to be considered, given as a list of indexes, list of atom names, 
or a selection following the syntax of `PDBTools`, or vector of `PDBTools.Atom`s, or a `PDBTools.Residue`

# Examples

In the following example we compute the coordination number of the atoms of residue 50 (of the solute) with the solvent atoms of TMAO,
as a function of the distance. Finally, we show the average number of TMAO molecules within 5 Angstroms of residue 50. 
The `findlast(<(5), R.d)` part of the code below returns the index of the last element of the `R.d` array that is smaller than 5 Angstroms.

## Precomputing the group contributions Using the `contrib` function

```julia
using ComplexMixtures, PDBTools
pdb = readPDB("test/data/NAMD/structure.pdb");
R = load("test/data/NAMD/protein_tmao.json");
solute = Selection(PDBTools.select(pdb, "protein"), nmols=1);
residue50 = PDBTools.select(pdb, "residue 50");
# Compute the group contributions to the MDDF
residue50_contribution = contrib(solute, R.solute_atom, residue50);
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
solute = Selection(PDBTools.select(pdb, "protein"), nmols=1);
residue50 = PDBTools.select(pdb, "residue 50");
# Compute the coordination number
residue50_coordination = coordination_number(solute, R.solute_atom, R, group)
# Output the average number of TMAO molecules within 5 Angstroms of residue 50
residue50_coordination[findlast(<(5), R.d)]
```

"""
function coordination_number(s::Selection, atom_contributions::Matrix{Float64}, R::Result, group)
    # Extract the group contributions to the MDDF
    group_contrib = contrib(s, atom_contributions, group)
    # Compute the coordination number
    return coordination_number(R, group_contrib)
end
function coordination_number(R::Result, group_contributions::Vector{Float64})
    cn = cumsum(group_contributions[i] * R.md_count_random[i] for i in eachindex(group_contributions)) 
    return cn
end


@testitem "coordination_number" begin
    using ComplexMixtures, PDBTools
    import ComplexMixtures.Testing: test_dir

    pdb = readPDB("$test_dir/data/NAMD/structure.pdb") 
    R = load("$test_dir/data/NAMD/protein_tmao.json")
    solute = Selection(PDBTools.select(pdb, "protein"), nmols=1)
    cn = coordination_number(solute, R.solute_atom, R, PDBTools.select(pdb, "residue 50"))
    @test cn[findlast(<(5), R.d)] ≈ 0.24999999999999997 atol=1e-10

    pdb = readPDB("$test_dir/data/NAMD/Protein_in_Glycerol/system.pdb")
    R = load("$test_dir/data/NAMD/Protein_in_Glycerol/protein_water.json") 
    group = PDBTools.select(pdb, "protein")
    solute = Selection(group, nmols=1)
    cn = coordination_number(solute, R.solute_atom, R, group)
    @test maximum(R.sum_md_count .- cn) < 1e-10
end