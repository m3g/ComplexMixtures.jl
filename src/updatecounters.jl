"""
    setbin(d,step)

$(INTERNAL)

Function that sets to which histogram bin a data point pertains simple, but important to keep consistency over all calls.

"""
setbin(d, step) = trunc(Int, d / step) + 1

"""
    itype(iatom::Int, natomspermol::Int)

$(INTERNAL)

Given the index of the atom in the vector of coordinates of the solute or the solvent,
returns the type of the atom, that is, the index of this atom within the molecule
(goes from 1 to natomspermol)

"""
function itype(iatom::Int, natomspermol::Int)
    itype = iatom % natomspermol
    return ifelse(itype == 0, natomspermol, itype)
end

# Calling using the structures of Solute and Solvent, to clear up the code above
itype(iatom::Int, s::Selection) = itype(iatom, s.natomspermol)

"""
    inbulk(d, R::Result)

$(INTERNAL)

Function that returns if a distance is in the bulk region or not, according to the options.

"""
inbulk(d, R::Result) = R.options.usecutoff ? (d >= R.dbulk && d < R.cutoff) : (d >= R.dbulk)

"""
    updatecounters!

$(INTERNAL)

Function that updates the counters in `R` and returns `n_dmin_in_bulk` given
the output of `cutoffdistances`.

If the solute and solvent selections are provided, 
update md_count, rdf_count and the atom-specific counters

returns: `n_dmin_in_bulk`: number of molecules with all the atoms in the bulk
         `n_dref_in_bulk`: number of molecules with the reference atom in the bulk

"""
function updatecounters!(R::Result, solute::Selection, solvent::Selection, system)

    # list of minimum distances
    list = system.list

    # Update the reference atom counter
    n_dref_in_bulk = 0
    for i = 1:solvent.nmols
        if list[i].ref_atom_within_cutoff
            ibin = setbin(list[i].d_ref_atom, R.options.binstep)
            R.rdf_count[ibin] += 1
        end
        n_dref_in_bulk += inbulk(list[i].d_ref_atom, R)
    end

    # Sort the vectors such that the elements with distances 
    # smaller than the cutoff are at the begining, this is used to random 
    # sample the bulk molecules afterwards
    partialsort_cutoff!(dmin_mol, R.cutoff, by = x -> x.d)

    # Add distances to the counters
    n_solvent_in_domain = 0
    i = 1
    while i <= solvent.nmols && dmin_mol[i].d < R.cutoff
        ibin = setbin(dmin_mol[i].d, R.options.binstep)
        R.md_count[ibin] += 1
        R.solute_atom[ibin, itype(dmin_mol[i].iat, solute)] += 1
        R.solvent_atom[ibin, itype(dmin_mol[i].jat, solvent)] += 1
        if !inbulk(dmin_mol[i].d, R)
            n_solvent_in_domain += 1
        end
        i = i + 1
    end
    if R.options.usecutoff
        n_dmin_in_bulk = (i - 1) - n_solvent_in_domain
    else
        n_dmin_in_bulk = solvent.nmols - (i - 1)
    end

    return n_dmin_in_bulk, n_dref_in_bulk
end

#
# If the rdf_count_random_frame is provided, update the counters associated to the random distribution
#
function updatecounters!(
    R::Result,
    rdf_count_random_frame::Vector{Float64},
    md_count_random_frame::Vector{Float64},
    solvent::Selection,
    system
)

    list = system.list
    # Update the reference atom counter
    for i = 1:solvent.nmols
        if list[i].ref_atom_within_cutoff 
            ibin = setbin(list[i].dref_mol, R.options.binstep)
            rdf_count_random_frame[ibin] += 1
        end
    end

    # Sort the vectors such that the elements with distances 
    # smaller than the cutoff are at the begining, this is used to random 
    # sample the bulk molecules afterwards
    partialsort_cutoff!(dmin_mol, R.cutoff, by = x -> x.d)

    # Add distances to the counters
    i = 1
    while i <= solvent.nmols && dmin_mol[i].d < R.cutoff
        ibin = setbin(list[i].dmin_mol, R.options.binstep)
        md_count_random_frame[ibin] += 1
        i = i + 1
    end

    return nothing
end
