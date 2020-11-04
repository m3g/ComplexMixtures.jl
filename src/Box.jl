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
  nc1 = max(1,trunc(Int64,sides[1]/(cutoff/lcell))) 
  nc2 = max(1,trunc(Int64,sides[2]/(cutoff/lcell))) 
  nc3 = max(1,trunc(Int64,sides[3]/(cutoff/lcell))) 
  nc = Vi3(nc1,nc2,nc3)
  l1 = sides[1]/nc1
  l2 = sides[2]/nc2
  l3 = sides[3]/nc3
  l = Vf3(l1,l2,l3) 
  Box(sides,nc,l,lcell)
end
