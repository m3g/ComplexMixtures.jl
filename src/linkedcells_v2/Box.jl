#
# Structure that contains some data required to compute the linked cells
#
struct Box
  cutoff :: Float64
  sides :: Vector{Float64}
  nc :: Vector{Int64}
end
# Must be initialized with a cutoff
Box(cutoff :: Float64) = Box( cutoff, zeros(3), zeros(Int64,3))
