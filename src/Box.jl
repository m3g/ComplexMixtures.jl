#
# Structure that contains some data required to compute the linked cells
#
struct Box
  sides :: Vector{Float64}
  nc :: Vector{Int64}
  l :: Vector{Float64}
  lcell :: Int64
end
# Must be initialized with the lcell parameters, which determines the fineness
# of the cell grid (cutoff/lcell)
Box(lcell) = Box(zeros(3), zeros(Int64,3), zeros(3), lcell)
