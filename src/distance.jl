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
  d = 0.
  for i in 1:3
    dx = (x[i]-y[i])%sides[i]
    if dx > sides[i]/2
      dx = dx - sides[i]
    elseif dx < -sides[i]/2
      dx = dx + sides[i]
    end
    d += dx^2
  end
  sqrt(d)
end
