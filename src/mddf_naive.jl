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

  # Computing all minimum-distances
  for iframe in 1:lastframe

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
    R.volume.total = R.volume.total + sides[1]*sides[2]*sides[3] 

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

      # cycle over solvent molecules
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
          # computing the number of solvent molecules in bulk, used for normalization
          # this normalization might be redudant if the solution is homogeneous (that is,
          # if the solute is not a single molecule at infinite dilution)
          if options.usecutoff 
            if dmin < cutoff
              n_jmol_inbulk = n_jmol_inbulk + 1
              jmol_inbulk[n_jmol_inbulk] = jmol
            end
          else
            n_jmol_inbulk = n_jmol_inbulk + 1
            jmol_inbulk[n_jmol_inbulk] = jmol
          end
        end
      end # solvent molecules 
      nbulk = nbulk + n_jmol_inbulk

      # Computing the random-solvent distribution to compute the random 
      # minimum-distance count
      for i in 1:nsamples
        # Choose randomly one molecule from the bulk
        jmol = rand(@view(jmol_inbulk[1:n_jmol_inbulk]))
        # Generate new random coordinates (translation and rotation) for this molecule
        jfmol = (jmol-1)*solvent.natomspermol + 1
        jlmol = jfmol + solvent.natomspermol - 1
        random_move!(jfmol,jlmol,x_solvent,sides,solute_center,x_solvent_random,moveaux)
        dmin, iatom, jatom, drefatom = minimumdistance(ifmol,ilmol,x_solute,1,solvent.natomspermol,x_solvent_random,options.irefatom)
        ibin = setbin(dmin,options.binstep)
        if ibin <= R.nbins
          R.count_random[ibin] += 1
        end
        # Use the position of the reference atom to compute the shell volume by Monte-Carlo integration
        ibin = setbin(drefatom,options.binstep)
        if ibin <= R.nbins
          R.volume.shell[ibin] += 1 
        end
      end

    end # solute molecules
  end # frames
  closetraj(trajectory)

  # Setup the distance vector
  for i in 1:R.nbins
    R.d[i] = shellradius(i,options.binstep)
  end

  #
  # Averaging
  #

  # Number of frames
  nframes = (lastframe - options.firstframe)/options.stride + 1 

  # Counters
  @. R.count = R.count / nframes
  @. R.count_random = R.count_random * solvent.nmols / (nsamples*nframes)
  @. R.solute_atom = R.solute_atom / nframes
  @. R.solvent_atom = R.solvent_atom / nframes

  # Volumes and Densities
  R.volume.total = R.volume.total / nframes
  R.density.solvent = solvent.nmols / R.volume.total
  R.density.solute = solute.nmols / R.volume.total

  R.volume.shell = R.volume.total * R.volume.shell / (nsamples*nframes)
  R.volume.domain = sum(R.volume.shell)
  R.volume.bulk = R.volume.total - R.volume.domain

  R.density.solvent_bulk = (nbulk/nframes) / R.volume.bulk 

  # Normalizing to compute distributions
  @. R.mddf = R.count / R.count_random
  for ibin in 1:R.nbins
    for i in 1:solute.natomspermol   
      R.solute_atom[i,ibin] = R.solute_atom[i,ibin] / R.count_random[ibin]
    end
    for j in 1:solvent.natomspermol
      R.solvent_atom[j,ibin] = R.solute_atom[j,ibin] / R.count_random[ibin]
    end
  end

  # Conversion factor for volumes (as KB integrals), from A^3 to cm^3/mol
  mole = 6.022140857e23
  convert = mole / 1.e24
  R.volume.total = convert * R.volume.total
  R.volume.shell = convert * R.volume.shell
  R.volume.domain = convert * R.volume.domain
  R.volume.bulk = convert * R.volume.bulk

  # Open output file and writes all information of this run

  #if print_files
  #  write_output_files(R)
  #end

  #if print_results
  #  print_results(R)
  #end

  return R

end

