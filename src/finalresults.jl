#
# Function that computes the final results of all the data computed by averaging
# according to the sampling of each type of data, and converts to common units
#
# Computes also the final distribution functions and KB integrals
#
# This function modified the values contained in the R data structure
#

struct Samples
  count :: Int64
  random :: Int64
end

function finalresults!(R :: Result, options :: Options, trajectory)
  
  # Conversion factor for volumes (as KB integrals), from A^3 to cm^3/mol
  mole = 6.022140857e23
  convert = mole / 1.e24

  # Setup the distance vector
  for i in 1:R.nbins
    R.d[i] = shellradius(i,options.binstep)
  end

  #
  # Averaging for the number of frames and number of solute molecules
  #
  s = Samples(R.nframes_read*trajectory.solute.nmols, 
              R.nframes_read*options.n_random_samples)

  # Counters
  @. R.md_count = R.md_count / s.count
  @. R.solute_atom = R.solute_atom / s.count
  @. R.solvent_atom = R.solvent_atom / s.count
  @. R.md_count_random = R.md_count_random / s.random
  @. R.rdf_count = R.rdf_count / s.count
  @. R.rdf_count_random = R.rdf_count_random / s.random

  # Volumes and Densities
  R.volume.total = R.volume.total / R.nframes_read
  R.density.solvent = R.density.solvent / R.nframes_read
  R.density.solute = R.density.solute / R.nframes_read

  R.volume.shell = R.volume.shell / R.nframes_read
  R.volume.domain = R.volume.domain / R.nframes_read
  R.volume.bulk = R.volume.bulk / R.nframes_read

  R.density.solvent_bulk = R.density.solvent_bulk / R.nframes_read

  # Fix the number of random samples using the bulk density
  if options.density_fix
    density_fix = R.density.solvent_bulk/R.density.solvent
    @. R.md_count_random = R.md_count_random * density_fix 
  end

  #
  # Computing the distribution functions and KB integrals, from the MDDF
  # and from the RDF
  #

  for ibin in 1:R.nbins

    # For the MDDF

    if R.md_count_random[ibin] > 0.
      R.mddf[ibin] = R.md_count[ibin] / R.md_count_random[ibin]
      for i in 1:trajectory.solute.natomspermol   
        R.solute_atom[i,ibin] = R.solute_atom[i,ibin] / R.md_count_random[ibin]
      end
      for j in 1:trajectory.solvent.natomspermol
        R.solvent_atom[j,ibin] = R.solvent_atom[j,ibin] / R.md_count_random[ibin]
      end
    end
    if ibin == 1
      R.sum_md_count[ibin] = R.md_count[ibin]
      R.sum_md_count_random[ibin] = R.md_count_random[ibin]
    else
      R.sum_md_count[ibin] = R.sum_md_count[ibin-1] + R.md_count[ibin]
      R.sum_md_count_random[ibin] = R.sum_md_count_random[ibin-1] + R.md_count_random[ibin]
    end
    R.kb[ibin] = convert*(1/R.density.solvent_bulk)*(R.sum_md_count[ibin] - R.sum_md_count_random[ibin])

    # For the RDF

    if R.rdf_count_random[ibin] > 0.
      R.rdf[ibin] = R.rdf_count[ibin] / (R.volume.shell[ibin]*R.density.solvent_bulk)
      #or
      #R.rdf[ibin] = R.rdf_count[ibin] / R.rdf_count_random[ibin] 
    end
    if ibin == 1
      R.sum_rdf_count[ibin] = R.rdf_count[ibin]
      R.sum_rdf_count_random[ibin] = R.rdf_count_random[ibin]
    else
      R.sum_rdf_count[ibin] = R.sum_rdf_count[ibin-1] + R.rdf_count[ibin]
      R.sum_rdf_count_random[ibin] = R.sum_rdf_count_random[ibin-1] + R.rdf_count_random[ibin]
    end
    R.kb_rdf[ibin] = convert*(1/R.density.solvent_bulk)*(R.sum_rdf_count[ibin] - R.sum_rdf_count_random[ibin])

  end

end
