#
# mddf_naive
#
# Computes the MDDF using the naive algorithm consisting of computing all distances
# useful for development purposes
#
# http://m3g.iqm.unicamp.br/
# http://github.com/m3g/MDDF
#

# Start with explicit optional input parameters

function mddf_naive(solute :: SoluteOrSolvent,
                    solvent :: SoluteOrSolvent,
                    trajectory, # The type of trajectory defines the functions used to read it
                    output_name :: String,
                   ;firstframe :: Int64 = 1,
                    lastframe :: Int64 = -1,
                    stride :: Int64 = 1,
                    periodic :: Bool = true, 
                    binstep :: Float64 = 0.2,
                    irefatom :: Int64 = 1,
                    dbulk :: Float64 = 10.,
                    nintegral :: Int64 = 10,
                    cutoff :: Float64 = -1.,
                    n_random_samples :: Int64 = -1,
                    print_files :: Bool = true,
                    print_results :: Bool = true
                   )

  input = InputDetails(firstframe = firstframe,
                       lastframe = lastrame,
                       stride = stride,
                       periodic = periodic,
                       binstep = binstep,
                       irefatom = irefatom,
                       dbulk = dbulk,
                       nintegral = nintegral,
                       cutoff = cutoff,
                       n_random_samples = n_random_samples,
                       print_files = print_files,
                       print_results = print_results,
                      )

  return mddf_naive(solute,solvent,trajectory,output_name,input)

end

# start providing the InputDetail data structure

function mddf_naive(solute :: SoluteOrSolvent,
                    solvent :: SoluteOrSolvent,
                    trajectory,
                    output_name :: String,
                    input :: InputDetails)
  
  # compute ibulk from dbulk (distance from which the solvent is considered bulk,
  # in the estimation of bulk density)

  if input.cutoff > 0. 
    usecutoff = true
    ibulk = round(Int64,dbulk/binstep) + 1
  else
    usecutoff = false
    ibulk = round(Int64,dbulk/binstep) + 1
    cutoff = input.dbulk
  end

  # Check if this is a single-solute or homogeneous solution situation: 
  if solute.nmols == 1
    single_solute = true
  else
    single_solute = false
  end

  # The number of random samples for numerical normalization
  if input.n_random_samples == -1
    nsamples = max(10*solvent.nmols,solute.nmols)
  end
  if ! single_solute
    nsamples = round(Int64,nsamples/solute.nmols)
  end

  # Initializing the structure that carries all data
  mddf = MDDF_Data(nbins,solute,solvent,output_name,input)

  # Last atom to be read from the trajectory files (actually this is relevant
  # only for the NamdDCD format, I think)
  natoms = length(solute)+length(solvent)

  # Vector to annotate the molecules that belong to the bulk solution
  jmol_inbulk = Vector{Int64}(undef,solvent.nmols)

  # Vector that will contain randomly generated solvent molecules
  x_solvent_random = Array{Float64}(undef,solvent.natomspermol,3)
  
  # Total number of bulk molecules
  nbulk = 0

  # Computing all minimum-distances

  open(trajectory)
  if trajectory.nframes < lastframe
    error(" The number of frames of the trajectory is smaller than user-defined lastframe ")
  end
  for iframe in 1:lastframe

    # reading coordinates of next frame
    nextframe!(trajectory,solute,solvent)
    if iframe < firstrame 
      continue
    end
    if iframe%stride != 0
      continue
    end

    # get pbc sides in this frame
    sides = getsides(trajectory,iframe)
    mddf.volume.total = mddf.volume.total + sides[1]*sides[2]*sides[3] 

    # Check if the cutoff is not too large considering the periodic cell size
    if cutoff > sides[1]/2. || cutoff > sides[2]/2. || cutoff > sides[3]/2.
      error("in MDDF: cutoff or dbulk > periodic_dimension/2 ")
    end

    # associate names for code clarity
    x_solute = trajectory.x_solute
    x_solvent = tajectory.x_solvent

    # computing the minimum distances, cycle over solute molecules
    for imol in 1:solute.nmols

      # first and last atoms of the current solute molecule
      ifmol = (imol-1)*solute.natomspermol + 1
      ilmol = ifmol + solute.natomspermol - 1

      # compute center of coordinates of solute molecule to wrap solvent coordinates around it
      solute_center = centerofcoordinates(@view(x_solute[ifmol:ilmol]))
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
        ibin = trunc(Int64,sqrt(dmin)/binstep)+1
        if ibin <= nbins
          mddf.count[ibin] += 1
          mddf.solute_atom[iatom,ibin] += 1 
          mddf.solvent_atom[jatom,ibin] += 1 
        else
          # computing the number of molecules in bulk, if at infinite dilution
# voltar: ver o que fazer quando não é diluição infita, deveria pular isto aqui,
# mas tenho que tomar cuidado com o cutoff
          if usecutoff 
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

      #
      # Estimating the shell volumes using Monte-Carlo integration
      #
      for i in 1:nsamples
        xrnd = -sizes/2 + rand(Float64,3)*sizes/2 + solute_center
        dmin = minimumdistance(xrnd,ifmol,ilmol,x_solute)
        ibin = trunc(Int64,sqrt(dmin)/binstep)+1
        if ibin <= nbins
          mddf.volume.shell[ibin] += 1
        end
      end

      # Computing the random-solvent distribution to compute the random 
      # minimum-distance count
      for i in 1:nsamples
        # Choose randomly one molecule from the bulk
        jmol = rand(@view(jmol_inbulk[1:n_jmol_inbulk]))
        # Generate new random coordinates (translation and rotation) for this molecule
        jfmol = (jmol-1)*solvent.natomspermol + 1
        jlmol = jfmol + solvent.natomspermol - 1
        random_move!(jfmol,jlmol,x_solvent,sizes,solute_center,x_solvent_random)
        dmin, iatom, jatom = minimumdistance(ifmol,ilmol,x_solute,1,solvent.natomspermol,x_solvent_random)
        ibin = trunc(Int64,sqrt(dmin)/binstep)+1
        if ibin <= nbins
          mddf.count_random[ibin] += 1
        end
      end

    end # solute molecules
  end # frames
  closetraj(trajectory)

  # Averaging
  density_fix = (bulkdensity*av_totalvolume)/nrsolvent_random
  mddf.count_random = density_fix * mddf_count_random

  framecount = round(Int64,(lastframe-firstframe+1)/stride)
  mddf.count = mddf.count / framecount
  mddf.count_random = (mddf.count_random / n_random_samples) / framecount
  mddf.solute_atom = mddf.solute_atom / framecount
  mddf.solvent_atom = mddf.solvent_atom / framecount

  # Normalizing to compute distributions
  @. mddf.mddf = mddf.count / mddf.count_random
  for ibin in 1:nbins
    for i in 1:solute.natomspermol   
      mddf.solute_atom[i,ibin] = mddf.solute_atom[i,ibin] / mddf.count_random[ibin]
    end
    for j in 1:solvent.natomspermol
      mddf.solvent_atom[j,ibin] = mddf.solute_atom[j,ibin] / mddf.count_random[ibin]
    end
  end

  mddf.density.solvent = solvent.nmols / mddf.volume.total
  mddf.density.solute = solute.nmols / mddf.volume.total

  mddf.volume.shell = mddf.volume.shell / n_random_samples / framecount 
  mddf.volume.total = mddf.volume.total / framecount
  mddf.volume.bulk = mddf.volume.total - sum(mddf.volume.shell)

  if single_solute
    mddf.density.bulk = nbulk / mddf.volume.bulk
  else
    mddf.density.bulk = mddf.density.solvent
  end

  # Conversion factor for volumes (as KB integrals), from A^3 to cm^3/mol
  mole = 6.022140857e23
  convert = mole / 1.e24

  # Open output file and writes all information of this run

  #if print_files
  #  write_output_files(mddf)
  #end

  #if print_results
  #  print_results(mddf)
  #end

  return mddf 

end

