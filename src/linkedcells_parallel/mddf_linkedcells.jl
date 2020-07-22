#     
# mddf_linkedcells
#
# Computes the MDDF using linked cells  
#
# http://m3g.iqm.unicamp.br/
# http://github.com/m3g/MDDF
#  

mutable struct FrameData

  trajectory # trajectory format
  volume_frame :: Volume
  rdf_count_random_frame :: Vector{Float64}
  sides :: Vector{Float64}
  box :: Box
  solute_center :: Vector{Float64}
  
  dc :: CutoffDistances
  dmin_mol :: Vector{DminMol}
  dref_mol :: Vector{Float64}
  x_solvent_random :: Array{Float64} 

  lc_solvent :: LinkedCells

  moveaux :: MoveAux
  nsamples :: Int64

end

function mddf_linkedcells(trajectory, options :: Options)  

  # Simplify code by assigning some shortened names
  solute = trajectory.solute
  solvent = trajectory.solvent
  x_solute = trajectory.x_solute
  x_solvent = trajectory.x_solvent

  # The number of random samples for numerical normalization
  nsamples = options.n_random_samples*solvent.nmols

  # Number of threads
  nthreads = Threads.nthreads()

  # Initializing the structure that carries the result per thread
  R = [ Result(trajectory,options) for i in 1:nthreads ]

  # Safe passing of frame counter to the threads
  tframe = zeros(Int64,nthreads)
  

  # Create data structures required for multithreading
  framedata = Vector{FrameData}(undef,nthreads)
  for ithread in 1:nthreads
    framedata[ithread] = FrameData(deepcopy(trajectory),                 # trajectory
                                   Volume(R[1].nbins),                   # volume_frame
                                   zeros(R[1].nbins),                    # rdf_count_random_frame
                                   zeros(3),                             # sides 
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

  # Skip initial frames if required
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
      @. framedata[ifree].trajectory.x_solute = trajectory.x_solute
      @. framedata[ifree].trajectory.x_solvent = trajectory.x_solvent
      @. framedata[ifree].sides = trajectory.sides
      iframe += 1
      iframe_read += 1 
      tframe[ifree] = iframe
      # Spawn the calculations for this frame
      t[ifree] = Threads.@spawn mddf_frame!(tframe[ifree],framedata[ifree],options,R[ifree])
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

  # Sum up the results of all threads into the data of thread one
  for ithread in 2:nthreads
    sum!(R[1],R[ithread])
  end

  # Setup the final data structure with final values averaged over the number of frames,
  # sampling, etc, and computes final distributions and integrals
  s = Samples(R[1].nframes_read*trajectory.solute.nmols,
              R[1].nframes_read*options.n_random_samples)
  finalresults!(R[1],options,trajectory,s)

  return R[1]

end

#
# Computes the MDDF for a single frame, modifies the data in the R (Result) structure
#

function mddf_frame!(iframe :: Int64, framedata :: FrameData, options :: Options, R :: Result)
  println(iframe)

  # Simplify code by assigning some shortened names
  trajectory = framedata.trajectory
  volume_frame = framedata.volume_frame
  rdf_count_random_frame = framedata.rdf_count_random_frame
  sides = framedata.sides
  box = framedata.box
  solute_center = framedata.solute_center
  dc = framedata.dc
  dmin_mol = framedata.dmin_mol
  dref_mol = framedata.dref_mol
  x_solvent_random = framedata.x_solvent_random
  lc_solvent = framedata.lc_solvent
  nsamples = framedata.nsamples
  moveaux = framedata.moveaux
  solute = trajectory.solute
  solvent = trajectory.solvent
  x_solute = trajectory.x_solute
  x_solvent = trajectory.x_solvent

  # Reset counters for this frame
  reset!(volume_frame)
  @. rdf_count_random_frame = 0.
  
  # get pbc sides in this frame
  sides = getsides(trajectory,iframe)

  volume_frame.total = sides[1]*sides[2]*sides[3]
  R.volume.total = R.volume.total + volume_frame.total
  
  R.density.solute = R.density.solute + (solute.nmols / volume_frame.total)
  R.density.solvent = R.density.solvent + (solvent.nmols / volume_frame.total)

  # Add the box side information to the box structure, in this frame
  @. box.sides = sides
  # Compute the number of cells in each dimension
  @. box.nc = max(1,trunc(Int64,box.sides/(options.cutoff/box.lcell)))
  @. box.l = box.sides/box.nc

  # Will wrap everthing relative to the reference atom of the first molecule
  # and move everything such that that center is in the origin. This is important
  # to simplify the computation of cell indexes, as the minimum coordinates are 
  # automatically -side/2 at each direction
  @. solute_center = x_solute[R.irefatom,1:3]
  wrap!(x_solute,sides,solute_center)
  center_to_origin!(x_solute,solute_center)
  wrap!(x_solvent,sides,solute_center)
  center_to_origin!(x_solvent,solute_center)

  # Initialize linked cells
  initcells!(x_solvent,box,lc_solvent)

  # Check if the cutoff is not too large considering the periodic cell size
  if options.cutoff > sides[1]/2. || options.cutoff > sides[2]/2. || options.cutoff > sides[3]/2.
    error("in MDDF: cutoff or dbulk > periodic_dimension/2 in frame: $iframe")
  end

  n_solvent_in_bulk = 0
  local n_solvent_in_bulk_last
  for isolute in 1:solute.nmols
    # We need to do this one solute molecule at a time to avoid exploding the memory
    # requirements
    x_this_solute = viewmol(isolute,x_solute,solute)

    # Compute all distances between solute and solvent atoms which are smaller than the 
    # cutoff (this is the most computationally expensive part), the distances are returned
    # in the dc structure
    cutoffdistances!(options.cutoff,x_this_solute,x_solvent,lc_solvent,box,dc)

    # For each solute molecule, update the counters (this is highly suboptimal, because
    # within updatecounters there are loops over solvent molecules, in such a way that
    # this will loop with cost nsolute*nsolvent. However, I cannot see an easy solution 
    # at this point with acceptable memory requirements
    n_solvent_in_bulk_last = updatecounters!(R,solute,solvent,dc,options,dmin_mol,dref_mol)
    n_solvent_in_bulk += n_solvent_in_bulk_last
  end

  #
  # Computing the random-solvent distribution to compute the random minimum-distance count
  #
  for i in 1:options.n_random_samples

    # generate random solvent box, and store it in x_solvent_random
    for j in 1:solvent.nmols
      # Choose randomly one molecule from the bulk, if there are actually bulk molecules
      if n_solvent_in_bulk_last != 0
        jmol = dmin_mol[rand(solvent.nmols-n_solvent_in_bulk_last+1:solvent.nmols)].jmol
      else
        jmol = rand(1:solvent.nmols)
      end
      # Indexes of this molecule in the x_solvent array
      x_ref = viewmol(jmol,x_solvent,solvent)
      # Indexes of the random molecule in random array
      x_rnd = viewmol(j,x_solvent_random,solvent)
      # Generate new random coordinates (translation and rotation) for this molecule
      random_move!(x_ref,R.irefatom,sides,x_rnd,moveaux)
    end

    # wrap random solvent coordinates to box, with the center at the origin
    wrap!(x_solvent_random,sides)

    # Initialize linked cells
    initcells!(x_solvent_random,box,lc_solvent)

    # Choose randomly one solute molecule to be the solute in this sample
    i_rand_mol = rand(1:solute.nmols)
    x_this_solute = viewmol(i_rand_mol,x_solute,solute)

    # Compute all distances between solute and solvent atoms which are smaller than the 
    # cutoff (this is the most computationally expensive part), the distances are returned
    # in the dc structure
    cutoffdistances!(options.cutoff,x_this_solute,x_solvent_random,lc_solvent,box,dc)

    # Update the counters of the random distribution
    updatecounters!(R.irefatom,R.md_count_random,rdf_count_random_frame,
                    solvent,dc,options,dmin_mol,dref_mol)


  end # random solvent sampling

  # Update counters with the data of this frame
  update_counters_frame!(R, rdf_count_random_frame, volume_frame, solute,
                         nsamples, n_solvent_in_bulk)

end

#
# Sum the counts of two Results structures, adding the result to the first structure
# R1 = R1 + R2
#

function sum!( R1 :: Result, R2 :: Result )

  @. R1.md_count += R2.md_count
  @. R1.md_count_random += R2.md_count_random

  @. R1.solute_atom += R2.solute_atom
  @. R1.solvent_atom += R2.solvent_atom

  @. R1.rdf_count += R2.rdf_count
  @. R1.rdf_count_random += R2.rdf_count_random

  sum!(R1.density,R2.density)
  sum!(R1.volume,R2.volume)

end

function sum!( D1 :: Density, D2 :: Density )
  D1.solute += D2.solute
  D1.solvent += D2.solvent
  D1.solvent_bulk += D2.solvent_bulk
end

function sum!( V1 :: Volume, V2 :: Volume )
  V1.total += V2.total
  V1.bulk += V2.bulk
  V1.domain += V2.domain
  @. V1.shell += V2.shell
end
