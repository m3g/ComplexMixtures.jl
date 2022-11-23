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

"""
    update_counters_frame!

$(INTERNAL)

Update some data accumulated for a frame.
"""
function update_counters_frame!(
    R::Result,
    rdf_count_random_frame::Vector{Float64},
    md_count_random_frame::Vector{Float64},
    volume_frame::Volume,
    solute::Selection,
    solvent::Selection,
    n_solvent_in_bulk::Float64,
)

    # Volume of each bin shell and of the solute domain
    @. volume_frame.shell =
        volume_frame.total *
        (rdf_count_random_frame / (R.options.n_random_samples * solvent.nmols))
    volume_frame.domain = sum(volume_frame.shell)

    # Volume and density of the bulk region in this frame
    if !R.options.usecutoff
        volume_frame.bulk = volume_frame.total - volume_frame.domain
    else
        ibulk = setbin(R.dbulk + 0.5 * R.options.binstep, R.options.binstep)
        for i = ibulk:R.nbins
            volume_frame.bulk = volume_frame.bulk + volume_frame.shell[i]
        end
    end

    bulk_density_at_frame = n_solvent_in_bulk / volume_frame.bulk

    # The number of "random" molecules generated must be that to fill the complete
    # volume with the density of the bulk. Now that we know the density of the bulk
    # and the volume of the frame, we can adjust that random count
    fillup_factor = bulk_density_at_frame * volume_frame.total / solvent.nmols
    @. R.rdf_count_random = R.rdf_count_random + rdf_count_random_frame * fillup_factor
    @. R.md_count_random = R.md_count_random + md_count_random_frame * fillup_factor

    # Update final counters, which have to be normalized by the number of frames 
    # at the end 
    @. R.volume.shell = R.volume.shell + volume_frame.shell
    R.volume.bulk = R.volume.bulk + volume_frame.bulk
    R.volume.domain = R.volume.domain + volume_frame.domain
    R.density.solvent_bulk = R.density.solvent_bulk + bulk_density_at_frame

    return R  
end