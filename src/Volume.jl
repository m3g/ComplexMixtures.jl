#
# Structures to contain the volume values obtained from calculations
#

mutable struct Volume
  shell :: Vector{Float64} 
  total :: Float64 
  bulk :: Float64 
  domain :: Float64
end

Volume(nbins :: Int64) = Volume( zeros(Float64,nbins), 0., 0., 0. )


