#
# Structures to contain the volume values obtained from calculations
#

using Parameters

@with_kw mutable struct Volume

  shell :: Vector{Float64} = 0.
  total :: Float64 = 0.
  bulk :: Float64 = 0.

end

