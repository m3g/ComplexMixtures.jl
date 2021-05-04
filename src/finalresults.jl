"""

```
finalresults!(R::Result, options::Options, trajectory::Trajectory, samples::Samples)
```

Function that computes the final results of all the data computed by averaging
according to the sampling of each type of data, and converts to common units.

Computes also the final distribution functions and KB integrals

This function modified the values contained in the R data structure


"""
function finalresults!(R::Result, options::Options, trajectory::Trajectory, samples::Samples)
  
  # Setup the distance vector
  for i in 1:R.nbins
    R.d[i] = shellradius(i,options.binstep)
  end

  # Counters
  @. R.md_count = R.md_count / (samples.md*R.nframes_read)
  @. R.solute_atom = R.solute_atom / (samples.md*R.nframes_read)
  @. R.solvent_atom = R.solvent_atom / (samples.md*R.nframes_read)
  @. R.md_count_random = R.md_count_random / (samples.random*R.nframes_read)
  @. R.rdf_count = R.rdf_count / (samples.md*R.nframes_read)
  @. R.rdf_count_random = R.rdf_count_random / (samples.random*R.nframes_read)

  # Volumes and Densities
  R.volume.total = R.volume.total / R.nframes_read
  R.density.solvent = R.density.solvent / R.nframes_read
  R.density.solute = R.density.solute / R.nframes_read

  R.volume.shell = R.volume.shell / R.nframes_read
  R.volume.domain = R.volume.domain / R.nframes_read
  R.volume.bulk = R.volume.bulk / R.nframes_read

  R.density.solvent_bulk = R.density.solvent_bulk / R.nframes_read

  #
  # Computing the distribution functions and KB integrals, from the MDDF
  # and from the RDF
  #

  for ibin in 1:R.nbins

    # For the MDDF

    if R.md_count_random[ibin] > 0.
      R.mddf[ibin] = R.md_count[ibin] / R.md_count_random[ibin]
      for i in 1:trajectory.solute.natomspermol   
        R.solute_atom[ibin,i] = R.solute_atom[ibin,i] / R.md_count_random[ibin]
      end
      for j in 1:trajectory.solvent.natomspermol
        R.solvent_atom[ibin,j] = R.solvent_atom[ibin,j] / R.md_count_random[ibin]
      end
    end
    if ibin == 1
      R.sum_md_count[ibin] = R.md_count[ibin]
      R.sum_md_count_random[ibin] = R.md_count_random[ibin]
    else
      R.sum_md_count[ibin] = R.sum_md_count[ibin-1] + R.md_count[ibin]
      R.sum_md_count_random[ibin] = R.sum_md_count_random[ibin-1] + R.md_count_random[ibin]
    end
    R.kb[ibin] = units.Angs3tocm3permol*
                 (1/R.density.solvent_bulk)*(R.sum_md_count[ibin] - R.sum_md_count_random[ibin])

    # For the RDF

    if R.rdf_count_random[ibin] > 0.
      R.rdf[ibin] = R.rdf_count[ibin] / R.rdf_count_random[ibin] 
    end
    if ibin == 1
      R.sum_rdf_count[ibin] = R.rdf_count[ibin]
      R.sum_rdf_count_random[ibin] = R.rdf_count_random[ibin]
    else
      R.sum_rdf_count[ibin] = R.sum_rdf_count[ibin-1] + R.rdf_count[ibin]
      R.sum_rdf_count_random[ibin] = R.sum_rdf_count_random[ibin-1] + R.rdf_count_random[ibin]
    end
    R.kb_rdf[ibin] = units.Angs3tocm3permol*
                     (1/R.density.solvent_bulk)*(R.sum_rdf_count[ibin] - R.sum_rdf_count_random[ibin])

  end

  nothing
end
