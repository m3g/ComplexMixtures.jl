# Computes the center of coordinates of a vector coor, and returns it in the vector
# cm, which must be provided previously allocated
function centerofcoordinates!(cm :: Vector{Float64}, coor :: Union{Array{Float64},SubArray{Float64}} )
  n = size(coor,1)
  @. cm = 0.
  for i in 1:n
    cm[1] = cm[1] + coor[i,1]
    cm[2] = cm[2] + coor[i,2]
    cm[3] = cm[3] + coor[i,3]
  end
  @. cm = cm / n
end

