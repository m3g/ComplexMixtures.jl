#
# mddf_naive
#
# Computes the MDDF using the naive algorithm consisting of computing all distances
# useful for development purposes
#
# http://m3g.iqm.unicamp.br/
# http://github.com/m3g/MDDF
#  

function mddf_naive(trajectory, options :: Options)  

  # Simplify code by assigning some shortened names
  solute = trajectory.solute
  solvent = trajectory.solvent
  x_solute = trajectory.x_solute
  x_solvent = trajectory.x_solvent
  
  # The number of random samples for numerical normalization
  nsamples = options.n_random_samples*solvent.nmols

  # Initializing the structure that carries all resuts
  R = Result(trajectory,options)

  # Vector to annotate the molecules that belong to the bulk solution
  jmol_inbulk = Vector{Int64}(undef,solvent.nmols)

  # Vector that will contain randomly generated solvent molecules
  x_solvent_random = Array{Float64}(undef,solvent.natomspermol,3)

  # Vector that wil contain the solute center of coordinates at each frame
  solute_center = Vector{Float64}(undef,3)

  # Auxiliary structure to random generation of solvent coordiantes
  moveaux = MoveAux()
  
  # Counter for the total number of bulk molecules
  nbulk = 0

  # Last frame to be considered
  if options.lastframe == -1 
    lastframe = trajectory.nframes
  else
    lastframe = options.lastframe
  end

  # Structure to organize counters for each frame only
  volume_frame = Volume(R.nbins)

  # Computing all minimum-distances
  for iframe in 1:lastframe

    # Reset counters for this frame
    reset!(volume_frame)

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

    # computing the minimum distances, cycle over solute molecules
    for imol in 1:solute.nmols

      # first and last atoms of the current solute molecule
      ifmol = (imol-1)*solute.natomspermol + 1
      ilmol = ifmol + solute.natomspermol - 1

      # compute center of coordinates of solute molecule to wrap solvent coordinates around it
      centerofcoordinates!(solute_center,@view(x_solute[ifmol:ilmol,:]))
      wrap!(x_solvent,sides,solute_center)

      # counter for the number of solvent molecules in bulk
      n_jmol_inbulk = 0

      #
      # cycle over solvent molecules to compute the MDDF count
      #
      for jmol in 1:solvent.nmols

        # first and last atoms of this solvent molecule
        jfmol = (jmol-1)*solvent.natomspermol + 1
        jlmol = jfmol + solvent.natomspermol - 1

        # Compute minimum distance 
        dmin, iatom, jatom = minimumdistance(ifmol,ilmol,x_solute,jfmol,jlmol,x_solvent)

        # Update histograms
        ibin = setbin(dmin,options.binstep)
        if ibin <= R.nbins
          R.count[ibin] += 1
          R.solute_atom[iatom,ibin] += 1 
          R.solvent_atom[jatom,ibin] += 1 
        else
          n_jmol_inbulk = n_jmol_inbulk + 1
          jmol_inbulk[n_jmol_inbulk] = jmol
        end

      end # solvent molecules 

      #
      # Computing the random-solvent distribution to compute the random minimum-distance count
      #
      for i in 1:nsamples
        # Choose randomly one molecule from the bulk
        jmol = rand(@view(jmol_inbulk[1:n_jmol_inbulk]))
        # Generate new random coordinates (translation and rotation) for this molecule
        jfmol = (jmol-1)*solvent.natomspermol + 1
        jlmol = jfmol + solvent.natomspermol - 1
        random_move!(jfmol,jlmol,x_solvent,sides,solute_center,x_solvent_random,moveaux)
        dmin, iatom, jatom, drefatom = minimumdistance(ifmol,ilmol,x_solute,
                                                       1,solvent.natomspermol,x_solvent_random,
                                                       options.irefatom)
        ibin = setbin(dmin,options.binstep)
        if ibin <= R.nbins
          R.count_random[ibin] += 1
        end
        # Use the position of the reference atom to compute the shell volume by Monte-Carlo integration
        ibin = setbin(drefatom,options.binstep)
        if ibin <= R.nbins
          volume_frame.shell[ibin] += 1
        else
          volume_frame.bulk += 1 
        end
      end
      @. volume_frame.shell = volume_frame.total * (volume_frame.shell / nsamples)
      volume_frame.bulk = volume_frame.total * (volume_frame.bulk / nsamples)

      @. R.volume.shell = R.volume.shell + volume_frame.shell
      R.volume.bulk = R.volume.bulk + volume_frame.bulk
      R.volume.domain = R.volume.domain + volume_frame.total - volume_frame.bulk
  
      R.density.solvent_bulk = R.density.solvent_bulk + n_jmol_inbulk / volume_frame.bulk

    end # solute molecules
  end # frames
  closetraj(trajectory)

  # Setup the distance vector
  for i in 1:R.nbins
    R.d[i] = shellradius(i,options.binstep)
  end

  #
  # Averaging for the number of frames
  #

  # Number of frames
  nframes = (lastframe - options.firstframe)/options.stride + 1 

  # Counters
  @. R.count = R.count / nframes
  @. R.solute_atom = R.solute_atom / nframes
  @. R.solvent_atom = R.solvent_atom / nframes
  @. R.count_random = R.count_random / (nframes*options.n_random_samples)

  # Volumes and Densities
  R.volume.total = R.volume.total / nframes
  R.density.solvent = R.density.solvent / nframes
  R.density.solute = R.density.solute / nframes

  R.volume.shell = R.volume.shell / nframes
  R.volume.domain = R.volume.domain / nframes
  R.volume.bulk = R.volume.bulk / nframes

  R.density.solvent_bulk = R.density.solvent_bulk / nframes

  # Fix the number of random samples using the bulk density
  if options.density_fix
    density_fix = R.density.solvent_bulk/R.density.solvent
    @. R.count_random = R.count_random * density_fix 
  end

  # Conversion factor for volumes (as KB integrals), from A^3 to cm^3/mol
  mole = 6.022140857e23
  convert = mole / 1.e24

  #
  # Computing the minimum distance distribution function normalized by
  # the random count or the density computed from the shell volumes
  #
  sum_count = 0.
  sum_count_random = 0.
  sum_nshell = 0.
  for ibin in 1:R.nbins
    if R.count_random[ibin] > 0.
      R.mddf[ibin] = R.count[ibin] / R.count_random[ibin]
      for i in 1:solute.natomspermol   
        R.solute_atom[i,ibin] = R.solute_atom[i,ibin] / R.count_random[ibin]
      end
      for j in 1:solvent.natomspermol
        R.solvent_atom[j,ibin] = R.solute_atom[j,ibin] / R.count_random[ibin]
      end
    end
    nshell = R.volume.shell[ibin]*R.density.solvent_bulk
    sum_nshell = sum_nshell + nshell
    if R.volume.shell[ibin] > 0.
      R.mddf_shell[ibin] = R.count[ibin] / nshell
    end
    sum_count = sum_count + R.count[ibin]
    sum_count_random = sum_count + R.count_random[ibin]
    # KB integral as a function of the distance
    R.kb[ibin] = convert*(1/R.density.solvent_bulk)*(sum_count - sum_count_random)
    R.kb_shell[ibin] = convert*(1/R.density.solvent_bulk)*(sum_count - sum_nshell)
  end

  return R

end

