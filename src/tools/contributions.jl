"""
    contributions(R::Result, group::Union{SoluteGroup,SolventGroup}; type = :mddf)

Returns the contributions of the atoms of the solute or solvent to the MDDF, coordiantion number or MD count.

"""
function contributions(
    R::Result, 
    group::Union{SoluteGroup,SolventGroup};
    type = :mddf, 
)
    if group isa SoluteGroup
        sol = R.solute
        group_count = R.solute_group_count
    elseif group isa SolventGroup
        sol = R.solvent
        group_count = R.solvent_group_count
    end

    sel_count = zeros(length(group_count[1]))

    # If the index of the groups was provided
    if !isnothing(group.group_index)
        igroup = group.group_index
        if igroup > length(group_count)
            throw(ArgumentError("Group $igroup not found in the group contribution count array."))
        end
        sel_count .= group_count[igroup]
    end

    # If the name of the group was provided
    if !isnothing(group.group_name)
        igroup = findfirst(==(group.group_name), sol.group_names)
        if isnothing(igroup) || igroup > length(group_count)
            throw(ArgumentError("Group $igroup not found in the group contribution count array."))
        end
        sel_count .= group_count[igroup]
    end

    # If, instead, atom indices or names were provided, sum over the contributions of the atoms.
    # This sum is different if the structure contanis one or more than one molecule.
    # If the structure has more than one molecule, than the indices are indices of 
    # the atoms *within* the molecule. 

    # Given atom inidices, sum over the contributions of the atoms
    if !isnothing(group.atom_indices) 
        if isempty(group.atom_indices)
            throw(ArgumentError("Group selection is empty"))
        end
        for i in group.atom_indices
            itype = findfirst(==(i), sol.indices)
            isnothing(itype) && throw(ArgumentError("Atom index $i not found in group data."))
            if sol.nmols > 1
                itype = atom_type(i, sol.natomspermol)
            end
            sel_count .+= group_count[itype]
        end
    end

    # Given atom names, sum over the contributions of the atoms
    if !isnothing(group.atom_names) 
        if isempty(group.atom_names)
            throw(ArgumentError("Group selection is empty"))
        end
        for name in group.atom_names
            itype = findfirst(==(name), sol.names)
            isnothing(itype) && throw(ArgumentError("Atom name $name not found in group data."))
            if sol.nmols > 1
                itype = atom_type(i, sol.natomspermol)
            end
            sel_count .+= group_count[itype]
        end
    end

    # Convert to the desired type
    if type == :mddf
        for i in eachindex(R.md_count_random, sel_count)
            if R.md_count_random[i] == 0.0
                sel_count[i] = 0.0
            else
                sel_count[i] /= R.md_count_random[i]
            end
        end
    elseif type == :coordination_number
        sel_count .= cumsum(sel_count)
    elseif type == :md_count
        # do nothign, already md_count
    else
        throw(ArgumentError("type must be :mddf (default), :coordination_number or :md_count"))
    end

    return sel_count
end

function warning_nmols_types()
    println("""
        WARNING: There is more than one molecule in this selection.
                 Contributions are summed over all atoms of the same type.
    """)
end

@testitem "contributions" begin
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

    # fetch contribution by atom name 
    @test contributions(results, SoluteGroup("N")) == zeros(500)
    @test contributions(results, SoluteGroup("N")) == zeros(500)

    C1_contributions = contributions(solute, SoluteGroup("name C1"))
    @test length(C1_contributions) == 500

    H33_contributions = contributions(solute, SoluteGroup("name H33"))
    @test length(H33_contributions) == 500

    # solvent contributions fetching
    N_contributions = contributions(solvent, SolventGroup("name N"))
    @test length(N_contributions) == 500

    C1_contributions = contributions(solvent, SolventGroup("name C1"))
    @test length(C1_contributions) == 500

    H33_contributions = contributions(solvent, SolventGroup("name H33"))
    @test length(H33_contributions) == 500

end
