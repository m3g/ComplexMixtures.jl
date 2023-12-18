"""
    contributions(s::AtomSelection, atom_contributions::Matrix{Float64}, selection)

Extract the contribution of a given atom type selection from the solute or solvent atomic contributions to the MDDF.

`s` here is the solute or solvent selection (type `ComplexMixtures.AtomSelection`)
`atom_contributions` is the `R.solute_atom` or `R.solvent_atom` arrays of the `Result` structure,
and the last argument is the selection of atoms from the solute to be considered, given as a list 
of indexes, list of atom names, vector of `PDBTools.Atom`s, or a `PDBTools.Residue`. 

## Extended help

For selections of one molecule, the function has an additional keyword option `first_atom_is_ref` that is `false` by default. 
If set to `true`, the index first atom of the selection is considered as a reference atom. For example if a solute has 100 atoms,
but its first atom in the PDB file is number 901, the selection of indexes `[1, 2, 3]` will refer to atoms
with indexes `[901, 902, 903]`.

"""
function contributions(
    s::AtomSelection,
    atom_contributions::Matrix{Float64},
    indexes::Vector{Int};
    first_atom_is_ref = false,
)
    nbins = size(atom_contributions, 1)
    c = zeros(nbins)
    # If the selection is a single molecule, the indexes can be anything (as they are the numbers printed
    # in the PDB file)
    if s.nmols == 1
        # If the first atom is a reference atom, the indexes are shifted by the index of the first atom
        if first_atom_is_ref
            first_index = first(s.index) - 1
        else
            first_index = 0
        end
        for it in indexes
            ind = findfirst(isequal(first_index + it), s.index)
            if isnothing(ind)
                error("Index $it of input list not found in selection indexes list.")
            end
            c += @view(atom_contributions[:, ind])
        end
    # If more than one molecule, the index must correspond to an atom within one molecule
    else
        for it in indexes
            if it > s.natomspermol
                error(
                    "The index list contains atoms with indexes greater than the number of atoms of one molecule.",
                )
            end
            c += @view(atom_contributions[:, it])
        end
    end
    return c
end

#
# If a list of atom names is provided
#
function contributions(
    s::AtomSelection,
    atom_contributions::Matrix{Float64},
    names::Vector{String},
)
    indexes = Vector{Int}(undef, 0)
    for name in names
        index = findall(isequal(name), s.names)
        if length(index) == 0
            error(" Atom in input list is not part of solvent (or solute): $name")
        end
        append!(indexes, index)
    end
    return contributions(s, atom_contributions, indexes; first_atom_is_ref = true)
end

#
# If a list of atoms of PDBTools.Atom is provided
#
function contributions(
    s::AtomSelection,
    atom_contributions::Matrix{Float64},
    atoms::Vector{PDBTools.Atom};
    warning = true,
)
    (warning && s.nmols > 1) && warning_nmols_types()
    indexes = PDBTools.index.(atoms)
    # Check which types of atoms belong to this selection
    selected_types = which_types(s, indexes, warning = warning)
    return contributions(s, atom_contributions, selected_types)
end

#
# If a residue of type PDBTools.Residue is provided
#
function contributions(
    s::AtomSelection,
    atom_contributions::Matrix{Float64},
    residue::PDBTools.Residue;
    warning = true,
)
    (warning && s.nmols > 1) && warning_nmols_types()
    indexes = collect(residue.range)
    # Check which types of atoms belong to this selection
    selected_types = which_types(s, indexes, warning = warning)
    return contributions(s, atom_contributions, selected_types)
end

function warning_nmols_types()
    println("""
        WARNING: There is more than one molecule in this selection.
                 Contributions are summed over all atoms of the same type.
    """)
end

@testitem "solute position" begin
    using ComplexMixtures
    using PDBTools
    using ComplexMixtures.Testing

    dir = "$(Testing.data_dir)/PDB"
    atoms = readPDB("$dir/trajectory.pdb", "model 1")

    solute = AtomSelection(select(atoms, "resname TMAO and resnum 1"), nmols = 1)
    solvent = AtomSelection(
        select(atoms, "resname TMAO and resnum 2 or resname TMAO and resnum 3"),
        nmols = 2,
    )

    traj = Trajectory("$dir/trajectory.pdb", solute, solvent, format = "PDBTraj")
    results = mddf(traj)

    # solute contributions fetching
    N_contributions = contributions(solute, results.solute_atom, ["N"])
    @test length(N_contributions) == 500

    C1_contributions = contributions(solute, results.solute_atom, ["C1"])
    @test length(C1_contributions) == 500

    H33_contributions = contributions(solute, results.solute_atom, ["H33"])
    @test length(H33_contributions) == 500

    # solvent contributions fetching
    N_contributions = contributions(solvent, results.solvent_atom, ["N"])
    @test length(N_contributions) == 500

    C1_contributions = contributions(solvent, results.solvent_atom, ["C1"])
    @test length(C1_contributions) == 500

    H33_contributions = contributions(solvent, results.solvent_atom, ["H33"])
    @test length(H33_contributions) == 500

end
