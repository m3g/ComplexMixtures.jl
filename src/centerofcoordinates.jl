function centerofcoordinates(coor :: Array{Float64})
  n = length(coor)
  cm = zeros(3)
  for i in 1:n
    cm[1] = cm[1] + coor[i,1]
    cm[2] = cm[2] + coor[i,2]
    cm[3] = cm[3] + coor[i,3]
  end
  @. cm = cm / n
  return cm
end

