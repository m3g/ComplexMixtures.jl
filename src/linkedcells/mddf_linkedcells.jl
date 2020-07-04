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
  
  # Counter for the total number of bulk molecules
  nbulk = 0

  # Structure to organize counters for each frame only
  volume_frame = Volume(R.nbins)
  rdf_count_random_frame = zeros(R.nbins)

  # Initialize the linked cell structures
  lc_solute = LinkedCells(solute.natoms)
  lc_solvent = LinkedCells(solvent.natoms)
 
  # Structure that contains the cutoff, sides, and xmin, to organize the linked
  # cell calculations
  box = Box(options.cutoff)

  # Structure that will contain the temporary useful information of all the  
  # distances found to the be smaller than the cutoff, and the corresponding
  # atom indexes. The vectors of this structure might be resized during the calculations
  dc = CutoffDistances(solvent.natoms)

  # Vectors used to parse the minimum distance data
  dmin_mol = zeros(solvent.nmols)
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

    # Check if the cutoff is not too large considering the periodic cell size
    if options.cutoff > sides[1]/2. || options.cutoff > sides[2]/2. || options.cutoff > sides[3]/2.
      error("in MDDF: cutoff or dbulk > periodic_dimension/2 ")
    end

    # Add the box side information to the box structure, in this frame
    @. box.sides = sides
    # Compute the number of cells in each dimension
    @. box.nc = trunc(Int64,box.sides/box.cutoff) + 1

    # Compute all distances between solute and solvent atoms which are smaller than the 
    # cutoff (this is the most computationally expensive part), the distances are returned
    # in the dc structure
    cutoffdistances!(x_solute,x_solvent,lc_solute,lc_solvent,box,dc)

    # Let us assume that there is only one solute molecule for a while, and
    # annotate which are the minimum distances of each solvent molecule to this solute
    # molecule
    @. dmin_mol = +Inf
    @. dref_mol = +Inf
    for i in 1:dc.nd[1]
      jmol = solvent.imol[dc.jat[i]]
      if dc.d[i] < dmin_mol[jmol]
        dmin_mol[jmol] = dc.d[i]
      end
      if itype(dc.jat[i],solvent) == R.irefatom
        dref_mol[jmol] = dc.d[i] 
      end
    end

    # Add distances to the counters
    n_solvent_in_bulk = 0
    @. solvent_in_bulk = 0
    for i in 1:solvent.nmols
      if dmin_mol[i] <= options.dbulk
        ibin = setbin(dmin_mol[i],options.binstep)
        R.md_count[ibin] += 1
        #R.solute_atom[itype(dc.iat[i],solute),ibin] += 1
        #R.solvent_atom[itype(dc.jat[i],solvent),ibin] += 1
      else
        n_solvent_in_bulk += 1
        solvent_in_bulk[n_solvent_in_bulk] = i
      end
      if dref_mol[i] <= options.dbulk
        ibin = setbin(dref_mol[i],options.binstep)
        R.rdf_count[ibin] += 1
      end
    end

    #
    # Computing the random-solvent distribution to compute the random minimum-distance count
    #
    for i in 1:options.n_random_samples
      
      # Choose randomly one solute molecule to be the reference solute for this
      # random solvent box
      i_rand_mol = rand(1:solute.nmols)
      ifmol = (i_rand_mol-1)*solute.natomspermol+1
      ilmol = ifmol + solute.natomspermol - 1
      centerofcoordinates!(solute_center, ifmol, ilmol, x_solute)

      # generate random solvent box, and store it in x_solvent_random
      for j in 1:solvent.nmols
        # Choose randomly one molecule from the bulk
        jmol = solvent_in_bulk[rand(1:n_solvent_in_bulk)]
        # Indexes of this molecule in the x_solvent array
        jfmol = (jmol-1)*solvent.natomspermol + 1
        jlmol = jfmol + solvent.natomspermol - 1
        # Indexes of the random molecule in random array
        jfstore = (j-1)*solvent.natomspermol + 1
        jlstore = jfstore + solvent.natomspermol - 1
        # Generate new random coordinates (translation and rotation) for this molecule
        #xrand = @view(x_solvent_random[jfstore:jlstore,1:3])
        random_move!(jfmol,jlmol,x_solvent,R.irefatom,sides,solute_center,
                     #xrand,moveaux)
                     jfstore,jlstore,x_solvent_random,moveaux)
      end

      # Compute all distances between solute and solvent atoms which are smaller than the 
      # cutoff (this is the most computationally expensive part), the distances are returned
      # in the dc structure
      cutoffdistances!(ifmol,ilmol,x_solute,x_solvent_random,lc_solute,lc_solvent,box,dc)

      # Let us assume that there is only one solute molecule for a while, and
      # annotate which are the minimum distances of each solvent molecule to this solute
      # molecule
      @. dmin_mol = +Inf
      @. dref_mol = +Inf
      for i in 1:dc.nd[1]
        jmol = solvent.imol[dc.jat[i]]
        if dc.d[i] < dmin_mol[jmol]
          dmin_mol[jmol] = dc.d[i]
        end
        if itype(dc.jat[i],solvent) == R.irefatom
          dref_mol[jmol] = dc.d[i] 
        end
      end

      # Add distances to the counters
      for i in 1:solvent.nmols
        if dmin_mol[i] <= options.dbulk
          ibin = setbin(dmin_mol[i],options.binstep)
          R.md_count_random[ibin] += 1
        end
        if dref_mol[i] <= options.dbulk
          ibin = setbin(dref_mol[i],options.binstep)
          rdf_count_random_frame[ibin] += 1
        end
      end

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

