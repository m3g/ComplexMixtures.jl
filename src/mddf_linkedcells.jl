#     
# mddf_linkedcells
#
# Computes the MDDF using linked cells, serial version
#  

function mddf_linkedcells(trajectory :: Trajectory, options :: Options, 
                          samples :: Samples, mddf_compute!)  

  # Initializing the structure that carries all results
  R = Result(trajectory,options)

  # Data structure to be passed to mddf_frame
  framedata = FrameData(trajectory,R)

  # Print some information about this run
  title(R,trajectory.solute,trajectory.solvent)

  # Computing all minimum-distances
  progress = Progress(R.lastframe_read-options.firstframe+1,1)
  for iframe in 1:R.lastframe_read

    # reading coordinates of next frame
    nextframe!(trajectory)
    if iframe < options.firstframe 
      continue
    end
    if iframe%options.stride != 0
      continue
    end
    mddf_compute!(iframe,framedata,options,R)   

    next!(progress)
  end # frames
  closetraj(trajectory)

  # Setup the final data structure with final values averaged over the number of frames,
  # sampling, etc, and computes final distributions and integrals
  finalresults!(R,options,trajectory,samples)
  println(bars)

  return R

end

