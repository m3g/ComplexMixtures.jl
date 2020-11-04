#
# Structure that contains some data required to compute the linked cells
#
struct Box
  sides :: Vf3
  nc :: Vi3
  l :: Vf3
  lcell :: Int64
end

function Box(lcell :: Int64, sides :: T, cutoff :: Float64) where T <: Vf3
  # Compute the number of cells in each dimension
  nc = Vi3(max.(1,trunc.(Int64,sides/(cutoff/lcell))))
  l = Vf3(sides ./ nc)
  Box(sides,nc,l,lcell)
end
