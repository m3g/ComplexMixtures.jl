#
#    atom_type(iatom::Int, natomspermol::Int)
#
# Given the index of the atom in the vector of coordinates of the solute or the solvent,
# returns the type of the atom, that is, the index of this atom within the molecule
# (goes from 1 to natomspermol)
#
atom_type(iatom::Int, natomspermol::Int) = mod1(iatom, natomspermol)

@testitem "atom_type" begin
    using ComplexMixtures: atom_type
    @test atom_type(1, 3) == 1
    @test atom_type(4, 3) == 1
    @test atom_type(5, 3) == 2
    @test atom_type(6, 3) == 3
end

# Function that updates the MD counters of the groups (or atoms) of the solute
# and solvent
function update_group_count!(group_count, ibin, iatom, frame_weight, sol::AtomSelection)
    itype = atom_type(iatom, sol.natomspermol)
    if isnothing(sol.group_atom_indices)
        group_count[itype][ibin] += frame_weight
    else
        for (igroup, group) in enumerate(sol.group_atom_indices)
            if itype in group
                group_count[igroup][ibin] += frame_weight
            end
        end
    end
    return group_count
end

@testitem "update_group_count!" begin
    using ComplexMixtures: update_group_count!
    solute = AtomSelection(1:6, 3, nothing)
    solvent = AtomSelection(7:12, 3, nothing)
    R = Result(solute, solvent)
end

#
#    updatecounters!(R::Result, system::AbstractPeriodicSystem)
#
# Function that updates the minimum-distance counters in `R`
#
function updatecounters!(R::Result, system::AbstractPeriodicSystem, frame_weight::Float64)
    for md in system.list
        !md.within_cutoff && continue
        ibin = setbin(md.d, R.options.binstep)
        R.md_count[ibin] += frame_weight
        update_group_count!(R.solute_group_count, ibin, md.i, frame_weight, R.solute)
        update_group_count!(R.solvent_group_count, ibin, md.j, frame_weight, R.solvent)
        if md.ref_atom_within_cutoff
            ibin = setbin(md.d_ref_atom, R.options.binstep)
            R.rdf_count[ibin] += frame_weight
        end
    end
    return R
end
# Update counters for the ideal gas distributions
function updatecounters!(R::Result, system::AbstractPeriodicSystem, frame_weight::Float64, ::Val{:random})
    for md in system.list
        !md.within_cutoff && continue
        ibin = setbin(md.d, R.options.binstep)
        R.md_count_random[ibin] += frame_weight
        if md.ref_atom_within_cutoff
            ibin = setbin(md.d_ref_atom, R.options.binstep)
            R.rdf_count_random[ibin] += frame_weight
        end
    end
    return R
end
