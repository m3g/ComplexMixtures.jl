#
# Structure to contain a list of all the distances smaller than
# the cutoff
#
struct CutoffDistances

  # To store all distances smaller than the cutoff
  d :: Vector{Float64} # All distances smaller than the cutoff
  iat :: Vector{Int64} # atom of the solute
  jat :: Vector{Int64} # atom of the solvent
  nd :: Vector{Int64} # Number of small distances (vector size = 1), just to be mutable

end

# Function that zeroes all the values in this structure

function reset!( c :: CutoffDistances )
  @. c.d = 0.
  @. c.iat = 0
  @. c.jat = 0
  nd[1] = 0
end


