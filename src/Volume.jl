#
# Structures to contain the volume values obtained from calculations
#

@with_kw mutable struct Volume
  total :: Float64 
  bulk :: Float64 
  domain :: Float64
  shell :: Vector{Float64} 
end

Volume(nbins :: Int64) = Volume( 0., 0., 0.,  zeros(Float64,nbins) )

function reset!( v :: Volume ) 
  v.total = 0.
  v.bulk = 0.
  v.domain = 0.
  @. v.shell = 0.
end

function Base.show( io :: IO, v :: Volume ) 
  n = length(v.shell)
  println(" Mean total box volume: $(v.total) ")
  println(" Mean bulk volume: $(v.bulk) ")
  println(" Mean solute domain volume: $(v.domain) ")
  println(" Volumes of first, medium, and last solvation shells: $(v.shell[1]), $(v.shell[round(Int64,n/2)]), $(v.shell[n])")
end


