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
  jmol_in_bulk = Vector{Int64}(undef,solvent.nmols)

  # Vector that will contain randomly generated solvent molecules (one at a time, 
  # in this naive algoritm)
  x_solvent_random = Array{Float64}(undef,solvent.natomspermol,3)

  # Vector that wil contain the solute center of coordinates at each frame
  solute_center = Vector{Float64}(undef,3)

  # Auxiliary structure to random generation of solvent coordiantes
  moveaux = MoveAux(solvent.natomspermol)
  
  # Structure to organize counters for each frame only
  volume_frame = Volume(R.nbins)
  rdf_count_random_frame = zeros(R.nbins)

  # Computing all minimum-distances
  progress = Progress(R.nframes_read*solute.nmols,1)
  for iframe in 1:R.lastframe_read

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

    # Counter for the cumulative number of solvent molecules found to be in bulk
    # relative to each solute molecule
    n_solvent_in_bulk = 0

    # computing the minimum distances, cycle over solute molecules
    for imol in 1:solute.nmols
      next!(progress)

      # first and last atoms of the current solute molecule
      ifmol = (imol-1)*solute.natomspermol + 1
      ilmol = ifmol + solute.natomspermol - 1

      # compute center of coordinates of solute molecule to wrap solvent coordinates around it
      centerofcoordinates!(solute_center,ifmol,ilmol,x_solute)
      wrap!(x_solvent,sides,solute_center)

      # counter for the number of solvent molecules in bulk for this solute molecule
      n_jmol_in_bulk = 0

      #
      # cycle over solvent molecules to compute the MDDF count
      #
      for jmol in 1:solvent.nmols

        # first and last atoms of this solvent molecule
        jfmol = (jmol-1)*solvent.natomspermol + 1
        jlmol = jfmol + solvent.natomspermol - 1

        # Compute minimum distance 
        dmin, iatom, jatom, drefatom = minimumdistance(ifmol,ilmol,x_solute,
                                                       jfmol,jlmol,x_solvent,
                                                       R.irefatom)

        # Update histogram for computation of MDDF
        ibin = setbin(dmin,options.binstep)
        if ibin <= R.nbins
          R.md_count[ibin] += 1
          R.solute_atom[iatom,ibin] += 1 
          R.solvent_atom[jatom,ibin] += 1 
        else
          n_jmol_in_bulk += 1
          jmol_in_bulk[n_jmol_in_bulk] = jmol
        end

        # Update histogram for the computation of the RDF
        ibin = setbin(drefatom,options.binstep) 
        if ibin <= R.nbins
          R.rdf_count[ibin] += 1
        end

      end # solvent molecules 

      # Sum up the number of solvent molecules found in bulk for this solute to the total 
      n_solvent_in_bulk = n_solvent_in_bulk + n_jmol_in_bulk

      #
      # Computing the random-solvent distribution to compute the random minimum-distance count
      #
      for i in 1:nsamples
        # Choose randomly one molecule from the bulk
        jmol = jmol_in_bulk[rand(1:n_jmol_in_bulk)]
        # Generate new random coordinates (translation and rotation) for this molecule
        jfmol = (jmol-1)*solvent.natomspermol + 1
        jlmol = jfmol + solvent.natomspermol - 1
        random_move!(jfmol,jlmol,x_solvent,R.irefatom,sides,solute_center,x_solvent_random,moveaux)
        dmin, iatom, jatom, drefatom = minimumdistance(ifmol,ilmol,x_solute,
                                                       1,solvent.natomspermol,x_solvent_random,
                                                       R.irefatom)

        if dmin <= options.dbulk
          ibin = setbin(dmin,options.binstep)
          R.md_count_random[ibin] += 1
        end
        # Use the position of the reference atom to compute the shell volume by Monte-Carlo integration
        if drefatom <= options.dbulk
          ibin = setbin(drefatom,options.binstep)
          rdf_count_random_frame[ibin] += 1
        end
      end # random solvent sampling

    end # solute molecules

    @. R.rdf_count_random = R.rdf_count_random + rdf_count_random_frame
    @. volume_frame.shell = volume_frame.total * (rdf_count_random_frame/(nsamples*solute.nmols))
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

