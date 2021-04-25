#     
# mddf_linkedcells
#
# Computes the MDDF using linked cells, serial version
#  

function mddf_linkedcells(trajectory::Trajectory, options::Options, 
                          samples::Samples, RNG, mddf_compute!)  

  # Initializing the structure that carries all results
  R = Result(trajectory,options)

  # Data structure to be passed to mddf_frame
  framedata = FrameData(trajectory,R)

  # Print some information about this run
  options.silent || title(R,trajectory.solute,trajectory.solvent)

  # Computing all minimum-distances
  if ! options.silent 
    progress = Progress(R.nframes_read,1)
  end
  for iframe in 1:R.lastframe_read

    # reading coordinates of next frame
    nextframe!(trajectory)
    if iframe < options.firstframe 
      continue
    end
    if iframe%options.stride != 0
      continue
    end
    mddf_compute!(iframe,framedata,options,RNG,R)   

    options.silent || next!(progress)
  end # frames
  closetraj(trajectory)

  # Setup the final data structure with final values averaged over the number of frames,
  # sampling, etc, and computes final distributions and integrals
  finalresults!(R,options,trajectory,samples)
  options.silent || println(bars)

  return R

end

