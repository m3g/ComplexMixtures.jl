#
# Functions to compute Euclidean distances
#

@inline function distance(x :: AbstractVector{Float64}, y :: AbstractVector{Float64})  
  return sqrt((y[1]-x[1])^2 + (y[2]-x[2])^2 + (y[3]-x[3])^2)
end

#
# With periodic conditions 
#

@inline function distance( sides :: Vector{Float64}, 
                           x :: AbstractVector{Float64}, y :: AbstractVector{Float64}) 
  # Wrap x relative to y
  wrapone!(x,sides,y)
  return distance(x,y)
end

