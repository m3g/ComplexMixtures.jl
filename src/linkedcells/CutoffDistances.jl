#
# Structure to contain a list of all the distances smaller than
# the cutoff
#
struct CutoffDistances

  # To store all distances smaller than the cutoff
  nd :: Vector{Int64} # Number of small distances (vector size = 1), just to be mutable
  d :: Vector{Float64} # All distances smaller than the cutoff
  iat :: Vector{Int64} # atom of the solute
  jat :: Vector{Int64} # atom of the solvent
  imol :: Vector{Int64} # molecule of the solute
  jmol :: Vector{Int64} # molecule of the solvent

  # Size of the arrays
  maxdim :: Vector{Int64}

end

# Generator

CutoffDistances( natoms :: Int64) =
  CutoffDistances(zeros(Int64,1), # nd
                  zeros(Float64,natoms), # d
                  zeros(Int64,natoms), # iat
                  zeros(Int64,natoms), # jat
                  zeros(Int64,natoms), # imol
                  zeros(Int64,natoms), # jmol
                  [ natoms ]) # maxdim

# Function that zeroes all the values in this structure

function reset!( c :: CutoffDistances )
  c.nd[1] = 0
  @. c.d = 0.
  @. c.iat = 0
  @. c.jat = 0
  @. c.imol = 0
  @. c.jmol = 0
end


