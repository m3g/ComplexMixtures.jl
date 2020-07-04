#
# Functions to compute Euclidean distances
#

distance(x :: AbstractVector{Float64}, y :: AbstractVector{Float64}) = 
  sqrt((y[1]-x[1])^2 + (y[2]-x[2])^2 + (y[3]-x[3])^2)

distance(x :: AbstractArray{Float64}, y :: AbstractArray{Float64}, j :: Int64) = 
  sqrt((y[j,1]-x[1])^2 + (y[j,2]-x[2])^2 + (y[j,3]-x[3])^2)

distance(x :: AbstractArray{Float64}, y :: AbstractArray{Float64}, i :: Int64, j :: Int64) = 
  sqrt((y[j,1]-x[i,1])^2 + (y[j,2]-x[i,2])^2 + (y[j,3]-x[i,3])^2)

# With possible periodic conditions 

function distance( box :: Box, x :: AbstractArray{Float64}, y :: AbstractArray{Float64}, 
                               i :: Int64, j :: Int64 )
 
  dx = (x[i,1]-y[j,1])%box.sides[1]
  dy = (x[i,2]-y[j,2])%box.sides[2]
  dz = (x[i,3]-y[j,3])%box.sides[3]

  if dx > box.sides[1]/2
    dx = dx - box.sides[1]
  elseif dx < -box.sides[1]/2
    dx = dx + box.sides[1]
  end

  if dy > box.sides[2]/2
    dy = dy - box.sides[2]
  elseif dy < -box.sides[2]/2 
    dy = dy + box.sides[2]
  end

  if dz > box.sides[3]/2
    dz = dz - box.sides[3]
  elseif dz < -box.sides[3]/2
    dx = dx + box.sides[3]
  end

  return sqrt(dx^2+dy^2+dz^2)

end

# To use with subarrarys

function distance( box :: Box, x :: AbstractArray{Float64}, y :: AbstractArray{Float64})
 
  dx = (x[1]-y[1])%box.sides[1]
  dy = (x[2]-y[2])%box.sides[2]
  dz = (x[3]-y[3])%box.sides[3]

  if dx > box.sides[1]/2
    dx = dx - box.sides[1]
  elseif dx < -box.sides[1]/2
    dx = dx + box.sides[1]
  end

  if dy > box.sides[2]/2
    dy = dy - box.sides[2]
  elseif dy < -box.sides[2]/2 
    dy = dy + box.sides[2]
  end

  if dz > box.sides[3]/2
    dz = dz - box.sides[3]
  elseif dz < -box.sides[3]/2
    dx = dx + box.sides[3]
  end

  return sqrt(dx^2+dy^2+dz^2)

end
