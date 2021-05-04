"""

```
merge(r::Vector{Result})
```

This function merges the results of MDDF calculations obtained by running the same
analysis on multiple trajectories, or multiple parts of the same trajectory. It returns
a Result structure of the same type, with all the functions and counters representing averages
of the set provided weighted by the number of frames read in each Result set.

"""
function merge(r::Vector{Result})

  nr = length(r)
  nframes_read = r[1].nframes_read
  error = false
  for ir in 2:nr
    nframes_read += r[ir].nframes_read
    if r[ir].nbins != r[1].nbins
      println("ERROR: To merge Results, the number of bins of the histograms of both sets must be the same.")
    end
    if (r[ir].cutoff - r[1].cutoff) > 1.e-8
      println("ERROR: To merge Results, cutoff distance of the of the histograms of both sets must be the same.")
    end
  end
  if error
    error(" Incompatible set of results to merge. ")
  end
  
  # List of files and weights
  nfiles = 0
  for ir in 1:nr
    nfiles += length(r[ir].files)
  end
  files = Vector{String}(undef,nfiles) 
  weights = Vector{Float64}(undef,nfiles)

  # Final resuls
  R = Result(options=r[1].options,
             nbins=r[1].nbins,
             dbulk=r[1].dbulk,
             cutoff=r[1].cutoff,
             irefatom=r[1].irefatom,
             lastframe_read=r[nr].lastframe_read,
             nframes_read=nframes_read,
             autocorrelation=r[1].autocorrelation,
             solute = r[1].solute,
             solvent = r[1].solvent,
             files=files,
             weights=weights) 

  # Average results weighting the data considering the number of frames of each data set
  @. R.d = r[1].d
  ifile = 0
  for ir in 1:nr
 
    w = r[ir].nframes_read / nframes_read

    @. R.mddf += w*r[ir].mddf
    @. R.kb += w*r[ir].kb

    @. R.rdf += w*r[ir].rdf
    @. R.kb_rdf += w*r[ir].kb_rdf

    @. R.md_count += w*r[ir].md_count 
    @. R.md_count_random += w*r[ir].md_count_random

    @. R.sum_md_count += w*r[ir].sum_md_count
    @. R.sum_md_count_random += w*r[ir].sum_md_count_random

    @. R.solute_atom += w*r[ir].solute_atom
    @. R.solvent_atom += w*r[ir].solvent_atom

    @. R.rdf_count += w*r[ir].rdf_count
    @. R.rdf_count_random += w*r[ir].rdf_count_random

    @. R.sum_rdf_count += w*r[ir].sum_rdf_count
    @. R.sum_rdf_count_random += w*r[ir].sum_rdf_count_random

    R.density.solute += w*r[ir].density.solute
    R.density.solvent += w*r[ir].density.solvent
    R.density.solvent_bulk += w*r[ir].density.solvent_bulk
 
    R.volume.total += w*r[ir].volume.total
    R.volume.bulk += w*r[ir].volume.bulk
    R.volume.domain += w*r[ir].volume.domain
    R.volume.shell += w*r[ir].volume.shell

    for j in 1:length(r[ir].files)
      ifile += 1
      R.files[ifile] = r[ir].files[j]
      R.weights[ifile] = w*r[ir].weights[j]
    end

  end

  return R

end
