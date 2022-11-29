"""
    itype(iatom::Int, natomspermol::Int)

$(INTERNAL)

Given the index of the atom in the vector of coordinates of the solute or the solvent,
returns the type of the atom, that is, the index of this atom within the molecule
(goes from 1 to natomspermol)

"""
itype(iatom::Int, natomspermol::Int) = mod1(iatom, natomspermol)

# Calling using the structures of Solute and Solvent, to clear up the code above
itype(iatom::Int, s::Union{SolSummary,Selection}) = itype(iatom, s.natomspermol)

"""
    updatecounters!(R::Result, system::AbstractPeriodicSystem)

$(INTERNAL)

Function that updates the minimum-distance counters in `R`

"""
function updatecounters!(R::Result, system::AbstractPeriodicSystem)
    for md in system.list
        !md.within_cutoff && continue
        ibin = setbin(md.d, R.options.binstep)
        R.md_count[ibin] += 1
        R.solute_atom[ibin, itype(md.i, R.solute)] += 1
        R.solvent_atom[ibin, itype(md.j, R.solvent)] += 1
        if md.ref_atom_within_cutoff
            ibin = setbin(md.d_ref_atom, R.options.binstep)
            R.rdf_count[ibin] += 1
        end
    end
    return R
end
# Update counters for the ideal gas distributions
function updatecounters!(R::Result, system::AbstractPeriodicSystem, ::Val{:random})
    for md in system.list
        !md.within_cutoff && continue
        ibin = setbin(md.d, R.options.binstep)
        R.md_count_random[ibin] += 1
        if md.ref_atom_within_cutoff
            ibin = setbin(md.d_ref_atom, R.options.binstep)
            R.rdf_count_random[ibin] += 1
        end
    end
    return R
end