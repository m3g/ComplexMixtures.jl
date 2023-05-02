"""
    coordination_number(s::Selection, atom_contributions::Array{Float64}, R::Result, group)

Computes the coordination number of a given group of atoms from the solute or solvent atomic contributions to the MDDF.

`s` here is the solute or solvent selection (type `ComplexMixtures.Selection`)

`atom_contributions` is the `R.solute_atom` or `R.solvent_atom` arrays of the `Result` structure

`R` is the `Result` structure,

and the last argument is the selection of atoms from the solute to be considered, given as a list of indexes, list of atom names, 
or a selection following the syntax of `PDBTools`, or vector of `PDBTools.Atom`s, or a `PDBTools.Residue`

# Example

In the following example we compute the coordination number of the solute atoms of residue 50 with the solvent atoms of TMAO,
as a function of the distance. Finally, we show the average number of TMAO molecules within 5 Angstroms of residue 50. 
The `findlast(<(5), R.d)` part of the code below returns the index of the last element of the `R.d` array that is smaller than 5 Angstroms.

```julia-repl
julia> using ComplexMixtures, PDBTools

julia> pdb = readPDB("test/data/NAMD/structure.pdb");

julia> R = load("test/data/NAMD/protein_tmao.json");

julia> solute = Selection(PDBTools.select(pdb, "protein"), nmols=1);

julia> cn = coordination_number(solute, R.solute_atom, R, PDBTools.select(pdb, "residue 50"))
500-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0
 ⋮
 0.24999999999999997
 0.24999999999999997
 0.24999999999999997
 0.24999999999999997

julia> cn[findlast(<(5), R.d)]
0.24999999999999997
```

"""
function coordination_number(s::Selection, atom_contributions::Array{Float64}, R::Result, group)
    # Extract the group contributions to the MDDF
    group_contrib = contrib(s, atom_contributions, group)
    # Compute the coordination number
    cn = cumsum(group_contrib[i] * R.md_count_random[i] for i in eachindex(group_contrib)) 
    return cn
end

@testitem "coordination_number" begin
    using ComplexMixtures, PDBTools
    using ComplexMixtures.Testing
    pdb = readPDB(pdbfile) 
    R = load("$test_dir/data/NAMD/protein_tmao.json")
    solute = Selection(PDBTools.select(pdb, "protein"), nmols=1)
    cn = coordination_number(solute, R.solute_atom, R, PDBTools.select(pdb, "residue 50"))
    @test cn[findlast(<(5), R.d)] ≈ 0.24999999999999997 atol=1e-14
end