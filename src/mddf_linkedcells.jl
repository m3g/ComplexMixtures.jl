#     
# mddf_linkedcells
#
# Computes the MDDF using linked cells  
#  

# With default input options

mddf_linkedcells(trajectory) = mddf_linkedcells(trajectory,Options())

# With explicit options provided

function mddf_linkedcells(trajectory, options :: Options)  

  # Simplify code by assigning some shortened names
  solute = trajectory.solute
  solvent = trajectory.solvent
  x_solute = trajectory.x_solute
  x_solvent = trajectory.x_solvent

  # Initializing the structure that carries all results
  R = Result(trajectory,options)

  # Check if the solute is the same as the solvent, and if so, use the self
  # routines to compute the mddf and normalize the data accordingly
  if solute.index != solvent.index
    mddf_compute! = mddf_frame!
    nsamples = options.n_random_samples*solvent.nmols
    s = Samples(R.nframes_read*trajectory.solute.nmols,
                R.nframes_read*options.n_random_samples)
  else
    mddf_compute! = mddf_frame_self!
    nsamples = options.n_random_samples
    npairs = round(Int64,solvent.nmols*(solvent.nmols-1)/2)
    nfix = solvent.nmols^2/npairs
    s = Samples(R.nframes_read*(trajectory.solvent.nmols-1),
                R.nframes_read*options.n_random_samples*nfix)
  end


  # The number of random samples for numerical normalization
  nsamples = options.n_random_samples*solvent.nmols

  # Data structure to be passed to mddf_frame
  framedata = FrameData(trajectory,                                       # trajectory
                        Volume(R.nbins),                                  # volume_frame
                        zeros(R.nbins),                                   # rdf_count_random_frame
                        Box(options.lcell),                               # box 
                        zeros(3),                                         # solute_center
                        CutoffDistances(solvent.natoms),                  # dc
                        [ DminMol(+Inf,i,0,0) for i in 1:solvent.nmols ], # dmin_mol
                        zeros(solvent.nmols),                             # dref_mol
                        similar(x_solvent),                               # x_solvent_random
                        LinkedCells(solvent.natoms),                      # lc_solvent
                        MoveAux(solvent.natomspermol),                    # moveaux
                        nsamples)                                         # nsamples        

  # Print some information about this run
  title(R,solute,solvent)

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
  finalresults!(R,options,trajectory,s)
  println(bars)

  return R

end

