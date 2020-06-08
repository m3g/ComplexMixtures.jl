#
# Structures to contain the results of the MDDF calculation
#

struct Result

  nbins :: Int64
  dmax :: Float64
  d :: Vector{Float64}
  mddf :: Vector{Float64}
  kb :: Vector{Float64}

  solute_atom :: Array{Float64}
  solvent_atom :: Array{Float64}
  
  count :: Vector{Float64}
  count_random :: Vector{Float64}
  
  density :: Density
  volume :: Volume

  file :: OutputFiles

  options :: Options

end

#
# Initialize the data structure that is returned from the computation, and checks some
# input parameters for consistency
#
              :
function Result( trajectory, options :: Options ) 

  # Names of auxiliary output files
  atom_contrib_solvent = FileOperations.remove_extension(options.output)*
                                 "-ATOM_CONTRIB_SOLVENT."*
                            FileOperations.file_extension(options.output)
  atom_contrib_solute = FileOperations.remove_extension(options.output)*
                            "-ATOM_CONTRIB_SOLUTE."*
                            FileOperations.file_extension(options.output)

  # Check for simple input errors
  if options.stride < 1
    error("in MDDF options: stride cannot be less than 1. ")
  end
  if options.lastframe > 0 && options.lastframe < options.firstframe
    error("in MDDF options: lastframe must be greater or equal to firstframe. ")
  end
  if options.lastframe > trajectory.nframes
    error("in MDDF options: lastframe is greater than trajectory.nframes. ")
  end

  # Check for problems in dbulk and cutoff definitions
  dmax = options.dbulk
  if options.usecutoff
    if options.dbulk >= options.cutoff 
      error(" in MDDF options: The bulk volume is zero (dbulk must be smaller than cutoff). ")
    end
    if (options.cutoff/options.binstep)%1 > 1.e-5
      error("in MDDF options: cutoff must be a multiple of binstep.")
    end
    if ((options.cutoff-options.dbulk)/options.binstep)%1 > 1.e-5
      error("in MDDF options: (cutoff-dbulk) must be a multiple of binstep. ")
    end
    dmax = options.cutoff
  end
  if (options.dbulk/options.binstep)%1 > 1.e-5
    error("in MDDF options: dbulk must be a multiple of binstep.")
  end
  nbins = round(Int64,dmax/options.binstep)

  if options.irefatom > trajectory.solvent.natoms 
    error("in MDDF options: Reference atom index", options.irefatom, " is greater than number of "*
          "                 atoms of the solvent molecule. ")
  end

  # Return data structure built up

  return Result( nbins, # number of bins of histogram
                 dmax, # maximum distance to be considered (cutoff or dbulk)
                 zeros(Float64,nbins), # d - vector of distances
                 zeros(Float64,nbins), # Vector to store the actual mddf
                 zeros(Float64,nbins), # Vector to store the KB integral
                 zeros(Float64,trajectory.solute.natomspermol,nbins), # Array to store the solute atom contributions
                 zeros(Float64,trajectory.solvent.natomspermol,nbins), # Array to store the solvent atom contributions
                 zeros(Float64,nbins), # Solvent count at each distance
                 zeros(Float64,nbins), # Random solvent count at each distance
                 Density(), # mutable scalars for results
                 Volume(nbins), # mutable vector and scalars for results
                 OutputFiles( options.output, # name of main output file
                              atom_contrib_solvent, # name of solvent atom contribution file
                              atom_contrib_solute ), # name of solute atom contribution file,
                 options # all input options
               )      

end

