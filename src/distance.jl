#
# Functions to compute Euclidean distances
#

@inline distance(x :: AbstractVector{Float64}, y :: AbstractVector{Float64}) =
  sqrt((y[1]-x[1])^2 + (y[2]-x[2])^2 + (y[3]-x[3])^2)

#
# With periodic conditions 
#

@inline function distance( x :: AbstractVector{Float64}, y :: AbstractVector{Float64}, 
                           sides :: Vector{Float64} )
  # Wrap x relative to y
  xnew = wrapone(x,sides,y)
  return distance(xnew,y)

end

