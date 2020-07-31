#     
# mddf_linkedcells
#
# Computes the MDDF using linked cells  
#
# http://m3g.iqm.unicamp.br/
# http://github.com/m3g/MDDF
#  

# With default input options

mddf_linkedcells_parallel(trajectory) = mddf_linkedcells_parallel(trajectory,Options())

# With explicit Options provided

function mddf_linkedcells_parallel(trajectory, options :: Options)  

  # Simplify code by assigning some shortened names
  solute = trajectory.solute
  solvent = trajectory.solvent
  x_solute = trajectory.x_solute
  x_solvent = trajectory.x_solvent

  # Number of threads
  nthreads = Threads.nthreads()

  # Initializing the structure that carries the result per thread
  R = [ Result(trajectory,options) for i in 1:nthreads ]

  # Check if the solute is the same as the solvent, and if so, use the self
  # routines to compute the mddf and normalize the data accordingly
  if solute.index != solvent.index
    mddf_compute! = mddf_frame!
    nsamples = options.n_random_samples*solvent.nmols
    s = Samples(R[1].nframes_read*trajectory.solute.nmols,
                R[1].nframes_read*options.n_random_samples)
  else
    mddf_compute! = mddf_frame_self!
    nsamples = options.n_random_samples
    npairs = round(Int64,solvent.nmols*(solvent.nmols-1)/2)
    nfix = solvent.nmols^2/npairs
    s = Samples(R[1].nframes_read*(trajectory.solvent.nmols-1),
                R[1].nframes_read*options.n_random_samples*nfix)
  end

  # Safe passing of frame counter to the threads
  tframe = zeros(Int64,nthreads)

  # Create data structures required for multithreading
  framedata = Vector{FrameData}(undef,nthreads)
  for ithread in 1:nthreads
    framedata[ithread] = FrameData(deepcopy(trajectory),                 # trajectory
                                   Volume(R[1].nbins),                   # volume_frame
                                   zeros(R[1].nbins),                    # rdf_count_random_frame
                                   Box(options.lcell),                   # box 
                                   zeros(3),                             # solute_center
                                   CutoffDistances(solvent.natoms),      # dc
                                   Vector{DminMol}(undef,solvent.nmols), # dmin_mol 
                                   zeros(solvent.nmols),                 # dref_mol
                                   similar(x_solvent),                   # x_solvent_random
                                   LinkedCells(solvent.natoms),          # lc_solvent
                                   MoveAux(solvent.natomspermol),        # moveaux
                                   nsamples)                             # nsamples
  end

  # Tasks
  t = Vector{Task}(undef,nthreads)
  free = ones(Bool,nthreads)

  # Skip initial frames if desired
  iframe = 0 # Counter for all frames of the input file
  while iframe < options.firstframe-1
    nextframe!(trajectory)
    iframe += 1
  end

  ndone = 0
  iframe_read = 0 # Counter for the frames that are actually being considered
  progress = Progress(R[1].nframes_read,1)
  while ndone < R[1].nframes_read

    # Launch for each free thread the computation of a frame
    while iframe_read < R[1].nframes_read && count(free) > 0 
      ifree = findfirst(x->x==true,free)
      # Skip frames that are not going to be considered 
      while (iframe+1)%options.stride != 0
        nextframe!(trajectory)
        iframe += 1
      end
      # Actually read the data of the frame that will be considered
      nextframe!(trajectory)
      iframe += 1
      iframe_read += 1 
      tframe[ifree] = iframe
      @. framedata[ifree].trajectory.x_solute = trajectory.x_solute
      @. framedata[ifree].trajectory.x_solvent = trajectory.x_solvent
      @. framedata[ifree].trajectory.sides = trajectory.sides
      # Spawn the calculations for this frame
      t[ifree] = Threads.@spawn mddf_compute!(tframe[ifree],framedata[ifree],options,R[ifree])
      free[ifree] = false
    end

    # Check thread status
    for ithread in 1:nthreads
      try 
        if istaskstarted(t[ithread])
          # If the task is done, add the results of this computation to
          # the final results
          if istaskdone(t[ithread])
            free[ithread] = true
            ndone += 1
            next!(progress)
          else
            free[ithread] = false
          end
        end
      catch 
        free[ithread] = true
      end
    end

    # Wait a little bit before checking all again
    sleep(options.sleep)

  end
  closetraj(trajectory)

  # Sum up the results of all threads into the data of thread one (R1<-R1+R2)
  for ithread in 2:nthreads
    sum!(R[1],R[ithread])
  end

  # Setup the final data structure with final values averaged over the number of frames,
  # sampling, etc, and computes final distributions and integrals
  finalresults!(R[1],options,trajectory,s)
  return R[1]

end

