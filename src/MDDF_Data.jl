#
# Structures to contain the results of the MDDF calculation
#

using Parameters

@with_kw mutable struct Density

  solute :: Float64 = 0.
  solvent :: Float64 = 0.
  solvent_bulk :: Float64 = 0.

end

@with_kw mutable struct Volume

  shell :: Vector{Float64} = 0.
  total :: Float64 = 0.
  bulk :: Float64 = 0.

end

@with_kw mutable struct OutputFiles

  output :: String
  solute_atoms :: String
  solvent_atoms :: String

end

@with_kw struct InputDetails

   firstframe :: Int64
   lastframe :: Int64
   stride :: Int64
   periodic :: Bool
   binstep :: Float64
   irefatom :: Int64
   dbulk :: Float64
   nintegral :: Int64
   cutoff :: Float64
   n_random_samples :: Int64
   print_files :: Bool
   print_results :: Bool

end

@with_kw struct MDDF_Data

  nbins :: Int64
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

  input :: InputDetails

end

# Initializer solute, solvent trajectory and output name, plus optional parameters

function MDDF_Data( solute :: Solute, 
                    solvent :: Solvent, 
                    trajectory,
                    output_name :: String, 
                    input :: InputDetails
                  ) 

  # Names of auxiliary output files
  atom_contrib_solvent = FileOperations.remove_extension(output_name)*
                                 "-ATOM_CONTRIB_SOLVENT."*
                            FileOperations.file_extension(output_name)
  atom_contrib_solute = FileOperations.remove_extension(output_name)*
                            "-ATOM_CONTRIB_SOLUTE."*
                            FileOperations.file_extension(output_name)

  # Check for simple input errors
  
  if input.stride < 1
    error("in MDDF input: stride cannot be less than 1. ")
  end
  if input.lastframe < input.firstframe && input.lastframe != 0
    error("in MDDF input: lastframe must be greater or equal to firstframe. ")
  end
  if input.dbulk - round(Int64,input.binstep*input.dbulk/input.binstep)*input.binstep > 1.e-5
    error("in MDDF input: dbulk must be a multiple of binstep.")
  end

  # compute ibulk from dbulk (distance from which the solvent is considered bulk,
  # in the estimation of bulk density)

  if input.cutoff > 0. 
    if input.dbulk >= input.cutoff 
      error(" in MDDF input: The bulk volume is zero (dbulk must be smaller than cutoff). ")
    end
    if (input.cutoff-input.dbulk)-round(Int64,(input.cutoff-input.dbulk)/input.binstep)*input.binstep > 1.e-5 
      error("in MDDF input: (cutoff-dbulk) must be a multiple of binstep. ")
    end
    nbins = round(Int64,input.cutoff/input.binstep)
  else
    println(" cutoff < 0: will use ntot-n(dbulk) as nbulk ")
    nbins = round(Int64,input.dbulk/input.binstep)
  end

  if input.irefatom > solvent.natoms 
    error("in MDDF input: Reference atom index", input.irefatom, " is greater than number of "*
          "               atoms of the solvent molecule. ")
  end

  # The number of random samples for numerical normalization
  if input.n_random_samples > 0 && input.n_random_samples < solute.nmols
    error("in MDDF input: n_random_samples must be greater of equal to the number of solute molecules. ")
  end

  # Return data structure built up

  return MDDF_Data( nbins = nbins, # number of bins of histogram
                    d = zeros(Float64,nbins),
                    mddf = Vector{Float64}(undef,nbins),
                    kb = Vector{Float64}(undef,nbins),
                    solute_atom = zeros(Float64,solute.natomspermol,nbins),
                    solvent_atom = zeros(Float64,solvent.natomspermol,nbins),
                    count = zeros(Float64,nbins),
                    count_random = zeros(Float64,nbins),
                    density = Density(), # mutable scalars for results
                    volume = Volume(), # mudtable scalars for results
                    file = OutputFiles( output_name, # name of main output file
                                        atom_contrib_solvent, # name of solvent atom contribution file
                                        atom_contrib_solute ), # name of solute atom contribution file,
                    input = input # all input options
                   )      

end


