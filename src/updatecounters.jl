#
# Function that updates the counters in R and returns n_dmin_in_bulk given
# the output of cutoffdistances
#
#

#
# If the solute and solvent selections are provided, 
# update md_count, rdf_count and the atom-specific counters
#
# returns: n_dmin_in_bulk: number of molecules with all the atoms in the bulk
#          n_dref_in_bulk: number of molecules with the reference atom in the bulk

function updatecounters!(R::Result, 
                         solute::Selection, solvent::Selection,
                         dc::CutoffDistances, 
                         dmin_mol::Vector{DminMol}, dref_mol::Vector{Float64})


  for i in 1:solvent.nmols
    dmin_mol[i].d = +Inf
    dref_mol[i] = +Inf
  end
  for i in 1:dc.nd[1]
    jmol = solvent.imol[dc.jat[i]]
    if dc.d[i] < dmin_mol[jmol].d
      dmin_mol[jmol].d = dc.d[i]
      dmin_mol[jmol].iat = dc.iat[i]
      dmin_mol[jmol].jat = dc.jat[i]
    end
    if itype(dc.jat[i],solvent) == R.irefatom
      if dc.d[i] < dref_mol[jmol]
        dref_mol[jmol] = dc.d[i] 
      end
    end
  end


  # Update the reference atom counter
  n_dref_in_bulk = 0
  for i in 1:solvent.nmols
    if dref_mol[i] < R.cutoff
      ibin = setbin(dref_mol[i],R.options.binstep)
      R.rdf_count[ibin] += 1
    end
    n_dref_in_bulk += inbulk(dref_mol[i],R)
  end

  # Sort the vectors such that the elements with distances 
  # smaller than the cutoff are at the begining, this is used to random 
  # sample the bulk molecules afterwards
  partialsort_cutoff!(dmin_mol,R.cutoff,by=x->x.d)

  # Add distances to the counters
  n_solvent_in_domain = 0
  i = 1
  while i <= solvent.nmols && dmin_mol[i].d < R.cutoff
    ibin = setbin(dmin_mol[i].d,R.options.binstep)
    R.md_count[ibin] += 1
    R.solute_atom[ibin,itype(dmin_mol[i].iat,solute)] += 1 
    R.solvent_atom[ibin,itype(dmin_mol[i].jat,solvent)] += 1 
    if ! inbulk(dmin_mol[i].d,R)
      n_solvent_in_domain += 1
    end
    i = i + 1
  end
  if R.options.usecutoff
    n_dmin_in_bulk = (i-1) - n_solvent_in_domain 
  else
    n_dmin_in_bulk = solvent.nmols - (i-1)
  end
  
  return n_dmin_in_bulk, n_dref_in_bulk
end

#
# If the rdf_count_random_frame is provided, update the
# counters associated to the random distribution
#
function updatecounters!(R::Result,
                         rdf_count_random_frame::Vector{Float64},
                         md_count_random_frame::Vector{Float64},
                         solvent::Selection, dc::CutoffDistances,
                         dmin_mol::Vector{DminMol}, dref_mol::Vector{Float64})

  for i in 1:solvent.nmols
    dmin_mol[i].d = +Inf
    dref_mol[i] = +Inf
  end
  for i in 1:dc.nd[1]
    jmol = solvent.imol[dc.jat[i]]
    if dc.d[i] < dmin_mol[jmol].d
      dmin_mol[jmol].d = dc.d[i]
    end
    if itype(dc.jat[i],solvent) == R.irefatom
      if dc.d[i] < dref_mol[jmol]
        dref_mol[jmol] = dc.d[i] 
      end
    end
  end

  # Update the reference atom counter
  for i in 1:solvent.nmols
    if dref_mol[i] < R.cutoff
      ibin = setbin(dref_mol[i],R.options.binstep)
      rdf_count_random_frame[ibin] += 1
    end
  end

  # Sort the vectors such that the elements with distances 
  # smaller than the cutoff are at the begining, this is used to random 
  # sample the bulk molecules afterwards
  partialsort_cutoff!(dmin_mol,R.cutoff,by=x->x.d)

  # Add distances to the counters
  i = 1
  while i <= solvent.nmols && dmin_mol[i].d < R.cutoff
    ibin = setbin(dmin_mol[i].d,R.options.binstep)
    md_count_random_frame[ibin] += 1
    i = i + 1
  end
  
  nothing
end

