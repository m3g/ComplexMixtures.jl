#     
# mddf_linkedcells
#
# Computes the MDDF using linked cells, in parallel
#

function mddf_linkedcells_parallel(trajectory :: Trajectory, options :: Options, 
                                   samples :: Samples, mddf_compute!)  

  # Simplify code by assigning some shortened names
  solute = trajectory.solute
  solvent = trajectory.solvent
  x_solute = trajectory.x_solute
  x_solvent = trajectory.x_solvent

  # Number of threads
  if options.nthreads < 0
    nspawn = Threads.nthreads() - 1
  else
    nspawn = options.nthreads - 1
  end
  if nspawn == 0
    error(" Parallel version must be executed only with more than 1 thread. ")
  end

  # Initializing the structure that carries the result per thread
  R0 = Result(trajectory,options)
  R = [ Result(trajectory,options,irefatom=R0.irefatom) for i in 1:nspawn ]

  # Safe passing of frame counter to the threads
  tframe = zeros(Int64,nspawn)

  # Create data structures required for multithreading
  framedata = Vector{FrameData}(undef,nspawn)
  for ispawn in 1:nspawn
    framedata[ispawn] = FrameData(deepcopy(trajectory),R[1])
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

  # Print some information about this run
  options.silent || title(R[1],solute,solvent,nspawn)

  ndone = 0
  iframe_read = 0 # Counter for the frames that are actually being considered
  if ! options.silent 
    progress = Progress(R[1].nframes_read,1)
  end
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
          options.silent || next!(progress)
          if options.GC && (Sys.free_memory() / Sys.total_memory() < options.GC_threshold)
            GC.gc() # why we need this anyway??? There should not be so much garbage.
          end
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
  finalresults!(R[1],options,trajectory,samples)
  options.silent || println(bars)

  return R[1]

end

