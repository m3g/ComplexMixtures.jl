#
# Structures to contain the results of the MDDF calculation
#

struct Result

  nbins :: Int64
  dmax :: Float64
  d :: Vector{Float64}

  # Data to compute the MDDF distribution and corresponding KB integral
  md_count :: Vector{Float64}
  md_count_random :: Vector{Float64}
  sum_md_count :: Vector{Float64}
  sum_md_count_random :: Vector{Float64}
  mddf :: Vector{Float64}
  kb :: Vector{Float64}

  # Atomic contributions to the MDDFs
  solute_atom :: Array{Float64}
  solvent_atom :: Array{Float64}
  
  # Data to compute a RDF and the KB integral from this count
  rdf_count :: Vector{Float64}
  rdf_count_random :: Vector{Float64}
  sum_rdf_count :: Vector{Float64}
  sum_rdf_count_random :: Vector{Float64}
  rdf :: Vector{Float64}
  kb_rdf :: Vector{Float64}

  # Overall densities and volumes
  density :: Density
  volume :: Volume

  # Options of the calculation
  options :: Options
  irefatom :: Int64
  lastframe_read :: Int64
  nframes_read :: Int64

  # File name(s) of the trajectories in this results 
  files :: Vector{String}
  weights :: Vector{Float64}

end

#
# Initialize the data structure that is returned from the computation, and checks some
# input parameters for consistency
#

function Result( trajectory, options :: Options ) 

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
  nbins = setbin(dmax,options.binstep)-1

  if options.irefatom > trajectory.solvent.natoms 
    error("in MDDF options: Reference atom index", options.irefatom, " is greater than number of "*
          "                 atoms of the solvent molecule. ")
  end

  # Set reference atom as the closest one to the center of coordinates of the molecule, as default
  if options.irefatom == -1
    nextframe!(trajectory)
    xfirst = @view(trajectory.x_solvent[1:trajectory.solvent.natomspermol,1:3])
    cm = centerofcoordinates(xfirst)
    dmin, one, irefatom = minimumdistance(cm,xfirst)
    firstframe(trajectory)
  else
    irefatom = options.irefatom
  end

  # Last frame to be considered
  if options.lastframe == -1
    lastframe_read = trajectory.nframes
  else
    lastframe_read = options.lastframe
  end
 
  # Actual number of frames that are read considering lastframe and stride
  nframes_read = round(Int64,(lastframe_read - options.firstframe + 1)/options.stride)

  # Return data structure built up

  return Result(options=options,
                nbins=nbins,
                dmax=dmax,
                irefatom=irefatom,
                lastframe_read=lastframe_read,
                nframes_read=nframes_read,
                solute_natomspermol=trajectory.solute.natomspermol,
                solvent_natomspermol=trajectory.solvent.natomspermol,
                files=[trajectory.filename],weights=[1.0])

end

#
# Generator with keyword parameters for the options that must be given
#

function Result(;options :: Options,
                nbins :: Int64,
                dmax :: Float64, 
                irefatom :: Int64, 
                lastframe_read :: Int64,
                nframes_read :: Int64,
                solute_natomspermol :: Int64,
                solvent_natomspermol :: Int64,
                files :: Vector{String},
                weights :: Vector{Float64})

  return Result( nbins, # number of bins of histogram
                 dmax, # maximum distance to be considered (cutoff or dbulk)
                 zeros(Float64,nbins), # d - vector of distances
                 zeros(Float64,nbins), # md_count
                 zeros(Float64,nbins), # md_count_random
                 zeros(Float64,nbins), # sum_md_count
                 zeros(Float64,nbins), # sum_md_count_random
                 zeros(Float64,nbins), # mddf
                 zeros(Float64,nbins), # kb
                 zeros(Float64,nbins,solute_natomspermol), # Array to store the solute atom contributions
                 zeros(Float64,nbins,solvent_natomspermol), # Array to store the solvent atom contributions
                 zeros(Float64,nbins), # rdf_count
                 zeros(Float64,nbins), # rdf_Count_random
                 zeros(Float64,nbins), # sum_rdf_count
                 zeros(Float64,nbins), # sum_rdf_count_random
                 zeros(Float64,nbins), # rdf
                 zeros(Float64,nbins), # kb_rdf
                 Density(), # mutable scalars for results
                 Volume(nbins), # mutable vector and scalars for results
                 options, # all input options
                 irefatom, # reference atom for RDF calculation
                 lastframe_read, # last frame read
                 nframes_read, # number of frames actually used for computing
                 files, # List of files that are considered in this results
                 weights # List of the weight of the data read from each file
               )      
end

#
# What to show at the REPL
#

function Base.show( io :: IO, R :: Result ) 

  println(" Trajectory files and weights: ")
  for i in 1:length(R.files) 
    println("   $(R.files[i]) - w = $(R.weights[i])")
  end
  println()

  ifar = trunc(Int64,R.nbins - 1.0/R.options.binstep)

  long_range_mean = mean( R.mddf[ifar:R.nbins] )
  long_range_std = std( R.mddf[ifar:R.nbins] )
  println(" Long range MDDF mean (expected 1.0): ", long_range_mean, " +/- ", long_range_std)

  long_range_mean = mean( R.rdf[ifar:R.nbins] )
  long_range_std = std( R.rdf[ifar:R.nbins] )
  println(" Long range RDF mean (expected 1.0): ", long_range_mean, " +/- ", long_range_std)

end




