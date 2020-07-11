#     
# mddf_linkedcells
#
# Computes the MDDF using linked cells  
#
# http://m3g.iqm.unicamp.br/
# http://github.com/m3g/MDDF
#  

function mddf_linkedcells(trajectory, options :: Options)  

  # Simplify code by assigning some shortened names
  solute = trajectory.solute
  solvent = trajectory.solvent
  x_solute = trajectory.x_solute
  x_solvent = trajectory.x_solvent

  # The number of random samples for numerical normalization
  nsamples = options.n_random_samples*solvent.nmols

  # Here, we have to generate the complete box of random solvent molecules at once,
  # to take advantadge of the linked cells, therefore we need an auxliary array
  # store that randomly generated solvent box
  x_solvent_random = similar(x_solvent)

  # Vector to annotate if the solvent molecule is a bulk molecule
  solvent_in_bulk = zeros(Int64,solvent.nmols)

  # Initializing the structure that carries all results
  R = Result(trajectory,options)

  # Vector that will contain the solute center of coordinates at each frame
  solute_center = Vector{Float64}(undef,3)

  # Vector that will store the reference-atom distances
  reference_distances = Vector{Float64}(undef,solvent.nmols)

  # Auxiliary structure to random generation of solvent coordinates
  moveaux = MoveAux(solvent.natomspermol)
  
  # Structure to organize counters for each frame only
  volume_frame = Volume(R.nbins)
  rdf_count_random_frame = zeros(R.nbins)

  # Initialize the linked cell structures
  lc_solvent = LinkedCells(solvent.natoms)
 
  # Structure that contains the sides, number of cells and cell length
  # in each dimension  to organize the linked cell calculations
  box = Box(options.lcell)

  # Structure that will contain the temporary useful information of all the  
  # distances found to the be smaller than the cutoff, and the corresponding
  # atom indexes. The vectors of this structure might be resized during the calculations
  dc = CutoffDistances(solvent.natoms)

  # Vectors used to parse the minimum distance data
  dmin_mol = Vector{DminMol}(undef,solvent.nmols)
  dref_mol = zeros(solvent.nmols)

  # Computing all minimum-distances
  progress = Progress(R.nframes_read,1)
  for iframe in 1:R.lastframe_read
    next!(progress)

    # Reset counters for this frame
    reset!(volume_frame)
    @. rdf_count_random_frame = 0.
  
    # reading coordinates of next frame
    nextframe!(trajectory)
    if iframe < options.firstframe 
      continue
    end
    if iframe%options.stride != 0
      continue
    end

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

    # Will wrap everthing relative to the center of coordinates of the solute atoms
    centerofcoordinates!(solute_center,x_solute)
    wrap!(x_solute,sides,solute_center)
    center_to_origin!(x_solute,solute_center)
    wrap!(x_solvent,sides,solute_center)
    center_to_origin!(x_solvent,solute_center)

    # Initialize linked cells
    initcells!(x_solvent,box,lc_solvent)

    # Check if the cutoff is not too large considering the periodic cell size
    if options.cutoff > sides[1]/2. || options.cutoff > sides[2]/2. || options.cutoff > sides[3]/2.
      error("in MDDF: cutoff or dbulk > periodic_dimension/2 ")
    end

    n_solvent_in_bulk = 0
    local n_solvent_in_bulk_last
    for isolute in 1:solute.nmols
      # We need to do this one solute molecule at a time to avoid exploding the memory
      # requirements
      ifmol = (isolute-1)*solute.natomspermol + 1
      ilmol = ifmol + solute.natomspermol - 1
      x_this_solute = @view(x_solute[ifmol:ilmol,1:3])

      # Compute all distances between solute and solvent atoms which are smaller than the 
      # cutoff (this is the most computationally expensive part), the distances are returned
      # in the dc structure
      cutoffdistances!(options.cutoff,x_this_solute,x_solvent,lc_solvent,box,dc)

      # For each solute molecule, update the counters (this is highly suboptimal, because
      # within updatecounters there are loops over solvent molecules, in such a way that
      # this will loop with cost nsolute*nsolvent. However, I cannot see an easy solution 
      # at this point with acceptable memory requirements
      n_solvent_in_bulk_last = updatecounters!(R.irefatom,R.md_count,R.rdf_count,
                                               solvent,dc,options,dmin_mol,dref_mol)
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
        jfmol = (jmol-1)*solvent.natomspermol + 1
        jlmol = jfmol + solvent.natomspermol - 1
        # Indexes of the random molecule in random array
        jfstore = (j-1)*solvent.natomspermol + 1
        jlstore = jfstore + solvent.natomspermol - 1
        # Generate new random coordinates (translation and rotation) for this molecule
        random_move!(@view(x_solvent[jfmol:jlmol,1:3]),R.irefatom,sides,
                     @view(x_solvent_random[jfstore:jlstore,1:3]),moveaux)
      end

      # wrap random solvent coordinates to box
      wrap!(x_solvent_random,sides,solute_center)
      center_to_origin!(x_solvent_random,solute_center)

      # Initialize linked cells
      initcells!(x_solvent_random,box,lc_solvent)

      # Choose randomly one solute molecule to be the solute in this sample
      i_rand_mol = rand(1:solute.nmols)
      ifmol = (i_rand_mol-1)*solute.natomspermol+1
      ilmol = ifmol + solute.natomspermol - 1
      x_this_solute = @view(x_solute[ifmol:ilmol,1:3])

      # Compute all distances between solute and solvent atoms which are smaller than the 
      # cutoff (this is the most computationally expensive part), the distances are returned
      # in the dc structure
      cutoffdistances!(options.cutoff,x_this_solute,x_solvent_random,lc_solvent,box,dc)

      # Update the counters and get the number of solvent molecules in bulk
      updatecounters!(R.irefatom,R.md_count_random,rdf_count_random_frame,
                      solvent,dc,options,dmin_mol,dref_mol)

    end # random solvent sampling

    @. R.rdf_count_random = R.rdf_count_random + rdf_count_random_frame
    @. volume_frame.shell = volume_frame.total * (rdf_count_random_frame/nsamples)
    volume_frame.domain = sum(volume_frame.shell)
    volume_frame.bulk = volume_frame.total - volume_frame.domain

    @. R.volume.shell = R.volume.shell + volume_frame.shell
    R.volume.bulk = R.volume.bulk + volume_frame.bulk
    R.volume.domain = R.volume.domain + volume_frame.domain
    R.density.solvent_bulk = R.density.solvent_bulk + (n_solvent_in_bulk/solute.nmols) / volume_frame.bulk

  end # frames
  closetraj(trajectory)

  # Setup the final data structure with final values averaged over the number of frames,
  # sampling, etc, and computes final distributions and integrals
  finalresults!(R,options,trajectory)

  return R

end

