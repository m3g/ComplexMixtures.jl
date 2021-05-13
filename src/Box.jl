"""

$(TYPEDEF)

Structure that contains some data required to compute the linked cells.

$(TYPEDFIELDS)

"""
struct Box{FloatVector,IntVector}
  sides::FloatVector
  nc::IntVector
  l::FloatVector
  lcell::Int
end

function Box(lcell::Int, sides::AbstractVector, cutoff::Float64)
  # Compute the number of cells in each dimension
  nc = SVector{3,Int}(max.(1,trunc.(Int,sides/(cutoff/lcell))))
  l = SVector{3,Float64}(sides ./ nc)
  Box(sides,nc,l,lcell)
end
