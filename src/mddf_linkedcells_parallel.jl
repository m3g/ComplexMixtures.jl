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

  # Initialize vector for random number generator
  init_random()

  # Simplify code by assigning some shortened names
  solute = trajectory.solute
  solvent = trajectory.solvent
  x_solute = trajectory.x_solute
  x_solvent = trajectory.x_solvent

  # Number of threads
  nspawn = Threads.nthreads()-1
  if nspawn == 0
    error(" Parallel version must be executed only with more than 1 thread. ")
  end


  # Initializing the structure that carries the result per thread
  R0 = Result(trajectory,options)
  R = [ Result(trajectory,options,irefatom=R0.irefatom) for i in 1:nspawn ]

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
  tframe = zeros(Int64,nspawn)

  # Create data structures required for multithreading
  framedata = Vector{FrameData}(undef,nspawn)
  for ispawn in 1:nspawn
    framedata[ispawn] = FrameData(deepcopy(trajectory),                 # trajectory
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
  t = Vector{Task}(undef,nspawn)
  free = ones(Bool,nspawn)

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
      t[ifree] = ThreadPools.@tspawnat ifree+1 mddf_compute!(tframe[ifree],framedata[ifree],options,R[ifree])
      free[ifree] = false
      #mddf_compute!(tframe[ifree],framedata[ifree],options,R[ifree])
    end

    # Wait a little bit before checking
    sleep(options.sleep)

    # Check thread status
    for ispawn in 1:nspawn
      if ! free[ispawn]
        if istaskfailed(t[ispawn])
          error(" Computation of MDDF failed in thread: $ispawn", fetch(t[ispawn]))
        end
        if istaskdone(t[ispawn])
          ndone += 1
          free[ispawn] = true
          next!(progress)
        end
      end
    end

  end
  closetraj(trajectory)

  # Sum up the results of all threads into the data of thread one (R1<-R1+R2)
  for ispawn in 2:nspawn
    sum!(R[1],R[ispawn])
  end

  # Setup the final data structure with final values averaged over the number of frames,
  # sampling, etc, and computes final distributions and integrals
  finalresults!(R[1],options,trajectory,s)
  return R[1]

end

