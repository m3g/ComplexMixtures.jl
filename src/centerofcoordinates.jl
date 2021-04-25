#
# Computes the center of coordinates of a vector coor
#
function centerofcoordinates(coor::AbstractVector{T}) where T
  cm = zeros(T)
  n = length(coor)
  for i in 1:n
    cm = cm + coor[i]
  end
  cm = cm / n
  return cm
end

