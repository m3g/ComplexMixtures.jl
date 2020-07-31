#
# This function merges the results of MDDF calculations obtained by running the same
# analysis on multiple trajectories, or multiple parts of the same trajectory. It returns
# a Result structure of the same type, with all the functions and counters representing averages
# of the set provided weighted by the number of frames read in each Result set.
#

function merge( results :: Vector{Results} )

  nr = length(results)
  nframes_read = r[1].nframes_read
  error = false
  for ir in 2:nr
    nframes_read += r.nframes_read
    if r[ir].nbins != r[1].nbins
      println("ERROR: To merge Results, the number of bins of the histograms of both sets must be the same.")
    end
    if (r[ir].dmax - r[1].dmax) > 1.d-8
      println("ERROR: To merge Results, the maximum distance of the of the histograms of both sets must be the same.")
    end
  end
  if error
    error(" Incompatible set of results to merge. ")
  end
  
  # Final resuls
  R = Result(options=r[1].options,
             nbins=r[1].nbins,
             dmax=r[1].dmax,
             irefatom=r[1].irefatom,
             lastframe_read=r[nr].lastframe_read,
             nframes_read=nframes_read) 

  # Average results weighting the data considering the number of frames of each data set
  
  for ir in 1:nr
 
    w = r[ir].nframes_read / nframes_read

    @. R.md_count = w*r[ir].md_count 
    @. R.md_count_random += w*r[ir].md_count_random

    @. R.solute_atom += w*r[ir].solute_atom
    @. R.solvent_atom += w*r[ir].solvent_atom

    @. R.rdf_count += w*r[ir].rdf_count
    @. R.rdf_count_random += w*r[ir].rdf_count_random

    R.density.solute += w*r[ir].density.solute
    R.density.solvent += w*r[ir].density.solvent
    R.density.solvent_bulk += w*r[ir].density.solvent_bulk
 
    R.volume.totall += w*r[ir].volume.total.
    R.volume.bulk += w*r[ir].volume.bulk
    R.volume.domain += w*r[ir].volume.domain
    R.volume.shell += w*r[ir].volume.shell

  end

  return R

end
