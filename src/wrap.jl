#
# Functions that wrap the coordinates (x{N,3} array) to obtain minimum images
# around a defined center
# 
# It modifies the coordinates input vector
#


function wrap!(x :: AbstractArray{Float64}, sides :: Vector{Float64}, center :: Vector{Float64})
  first = 1
  last = size(x,1)
  wrap!(first, last, x, sides, center)
end

function wrap!(first :: Int64, last :: Int64, x :: AbstractArray{Float64}, 
               sides :: Vector{Float64}, center :: Vector{Float64})
  for i in first:last
    wrapone!(i,x,sides,center)
  end
end

function wrapone!(i :: Int64, x :: AbstractArray{Float64}, sides :: Vector{Float64} ,center :: Vector{Float64})

  x[i,1] = x[i,1] - center[1]
  x[i,2] = x[i,2] - center[2]
  x[i,3] = x[i,3] - center[3]

  # Wrap to origin
  wrapone!(i,x,sides)

  x[i,1] = x[i,1] + center[1]
  x[i,2] = x[i,2] + center[2]
  x[i,3] = x[i,3] + center[3]

end

#
# If the center is not provided, wrap to origin
#


function wrap!(x:: AbstractArray{Float64}, sides :: Vector{Float64})
  first = 1
  last = size(x,1)
  wrap!(first,last,x,sides)
end

function wrap!(first :: Int64, last :: Int64, x :: AbstractArray{Float64}, sides :: Vector{Float64})
  for i in first:last
     wrapone!(i,x,sides)
  end
end

function wrapone!(i :: Int64, x :: AbstractArray{Float64}, sides :: Vector{Float64})

  x[i,1] = x[i,1]%sides[1]
  x[i,2] = x[i,2]%sides[2]
  x[i,3] = x[i,3]%sides[3]

  if x[i,1] > sides[1]/2  
    x[i,1] = x[i,1] - sides[1] 
  elseif x[i,1] < -sides[1]/2 
    x[i,1] = x[i,1] + sides[1]
  end

  if x[i,2] > sides[2]/2 
    x[i,2] = x[i,2] - sides[2] 
  elseif x[i,2] < -sides[2]/2 
    x[i,2] = x[i,2] + sides[2] 
  end

  if x[i,3] > sides[3]/2 
    x[i,3] = x[i,3] - sides[3]
  elseif x[i,3] < -sides[3]/2 
    x[i,3] = x[i,3] + sides[3]
  end

end










