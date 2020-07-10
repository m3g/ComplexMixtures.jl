#
# Function that updates the counters in R and returns n_solvent_in_bulk given
# the output of cutoffdistances
#

function updatecounters!(irefatom,md_count,rdf_count,solvent,dc,options,dmin_mol,dref_mol)

  for i in 1:solvent.nmols
    dmin_mol[i] = DminMol(+Inf,i)
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
      rdf_count[ibin] += 1
    end
  end

  # Sort the vectors such that the elements with distances 
  # smaller than the cutoff are at the begining, this is used to random 
  # sample the bulk molecules afterwards
  partialsort_cutoff!(dmin_mol,options.dbulk,by=x->x.d)

  # Add distances to the counters
  n_solvent_in_bulk = 0
  i = 1
  while i <= solvent.nmols && dmin_mol[i].d < options.dbulk 
    ibin = setbin(dmin_mol[i].d,options.binstep)
    md_count[ibin] += 1
    i = i + 1
  end
  n_solvent_in_bulk = solvent.nmols - i + 1

  return n_solvent_in_bulk

end
