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
  sort!(dmin_mol, by = mol -> mol.d)
  sort!(dref_mol)

  # Add distances to the counters
  n_solvent_in_bulk = 0
  i = 1
  while dmin_mol[i].d < options.dbulk 
    ibin = setbin(dmin_mol[i].d,options.binstep)
    md_count[ibin] += 1
    i = i + 1
  end
  n_solvent_in_bulk = solvent.nmols - i + 1
  i = 1
  while dref_mol[i] < options.dbulk
    ibin = setbin(dref_mol[i],options.binstep)
    rdf_count[ibin] += 1
    i = i + 1
  end

  return n_solvent_in_bulk

end
