#
# Function that updates the counters in R and returns n_solvent_in_bulk given
# the output of cutoffdistances
#

#
# If the R structure is provided, update md_count, rdf_count and the atom-specific
# counters
#

function updatecounters!(R :: Result, 
                         solute :: Selection, solvent :: Selection,
                         dc :: CutoffDistances, options :: Options,
                         dmin_mol :: Vector{DminMol}, dref_mol :: AbstractVector{Float64})

  for i in 1:solvent.nmols
    dmin_mol[i] = DminMol(+Inf,i,0,0)
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
  for i in 1:solvent.nmols
    if dref_mol[i] < R.cutoff
      ibin = setbin(dref_mol[i],options.binstep)
      R.rdf_count[ibin] += 1
    end
  end

  # Sort the vectors such that the elements with distances 
  # smaller than the cutoff are at the begining, this is used to random 
  # sample the bulk molecules afterwards
  partialsort_cutoff!(dmin_mol,R.cutoff,by=x->x.d)

  # Add distances to the counters
  n_solvent_in_bulk = 0
  i = 1
  while i <= solvent.nmols && dmin_mol[i].d < R.cutoff
    ibin = setbin(dmin_mol[i].d,options.binstep)
    R.md_count[ibin] += 1
    R.solute_atom[ibin,itype(dmin_mol[i].iat,solute)] += 1 
    R.solvent_atom[ibin,itype(dmin_mol[i].jat,solvent)] += 1 
    i = i + 1
    if inbulk(dmin_mol[i].d,R)
      n_solvent_in_bulk += 1
    end
  end
  return n_solvent_in_bulk

end

#
# If the md_count_random and rdf_count_random_frame are provided instead, update the
# counters associated to the random distribution
#

function updatecounters!(irefatom :: Int64, md_count_random :: AbstractVector{Float64},
                         rdf_count_random_frame :: AbstractVector{Float64},
                         solvent :: Selection, dc :: CutoffDistances,
                         options :: Options, 
                         dmin_mol :: Vector{DminMol}, dref_mol :: AbstractVector{Float64})

  for i in 1:solvent.nmols
    dmin_mol[i] = DminMol(+Inf,i,0,0)
    dref_mol[i] = +Inf
  end
  for i in 1:dc.nd[1]
    jmol = solvent.imol[dc.jat[i]]
    if dc.d[i] < dmin_mol[jmol].d
      dmin_mol[jmol].d = dc.d[i]
    end
    if itype(dc.jat[i],solvent) == irefatom
      if dc.d[i] < dref_mol[jmol]
        dref_mol[jmol] = dc.d[i] 
      end
    end
  end

  # Update the reference atom counter
  for i in 1:solvent.nmols
    if dref_mol[i] < options.dbulk
      ibin = setbin(dref_mol[i],options.binstep)
      rdf_count_random_frame[ibin] += 1
    end
  end

  # Sort the vectors such that the elements with distances 
  # smaller than the cutoff are at the begining, this is used to random 
  # sample the bulk molecules afterwards
  partialsort_cutoff!(dmin_mol,options.dbulk,by=x->x.d)

  # Add distances to the counters
  i = 1
  while i <= solvent.nmols && dmin_mol[i].d <= options.dbulk 
    ibin = setbin(dmin_mol[i].d,options.binstep)
    md_count_random[ibin] += 1
    i = i + 1
  end

end

