#
# Function that wraps the coordinates to obtain minimum images
# around a defined center
# 
# It modifies the coordinates input vector
#

function wrap!(coordinates :: Array{Float64}, sides :: Vector{Float64}, center :: Vector{Float64})
  for i in 1:length(coordinates)
    wrapone!(i,coordinates,sides,center)
  end
end

function wrapone!(i,x,sides,center)

  x[i,1] = x[i,1] - center[1]
  x[i,2] = x[i,2] - center[2]
  x[i,3] = x[i,3] - center[3]

  x[i,1] = x[i,1]%sides[1]
  x[i,2] = x[i,2]%sides[2]
  x[i,3] = x[i,3]%sides[3]

  if x[i,1] > sides[1]/2 ; x[i,1] = x[i,1] - sides[1] ; end
  if x[i,2] > sides[2]/2 ; x[i,2] = x[i,2] - sides[2] ; end
  if x[i,3] > sides[3]/2 ; x[i,3] = x[i,3] - sides[3] ; end

  if x[i,1] < -sides[1]/2 ; x[i,1] = x[i,1] + sides[1] ; end
  if x[i,2] < -sides[2]/2 ; x[i,2] = x[i,2] + sides[2] ; end
  if x[i,3] < -sides[3]/2 ; x[i,3] = x[i,3] + sides[3] ; end

  x[i,1] = x[i,1] + center[1]
  x[i,2] = x[i,2] + center[2]
  x[i,3] = x[i,3] + center[3]

end
