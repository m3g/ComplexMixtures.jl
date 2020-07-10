#
# Functions that wrap the coordinates (x{N,3} array) to obtain minimum images
# around a defined center
# 
# It modifies the coordinates input vector
#


function wrap!(x :: AbstractArray{Float64}, 
               sides :: Vector{Float64}, 
               center :: AbstractVector{Float64})
  for i in 1:size(x,1)
    wrapone!(@view(x[i,1:3]),sides,center)
  end
end

@inline function wrapone!(x :: AbstractVector{Float64}, 
                          sides :: Vector{Float64}, 
                          center :: AbstractVector{Float64})
  for i in 1:3
    x[i] = (x[i]-center[i])%sides[i]
    if x[i] > sides[i]/2  
      x[i] = x[i] - sides[i] 
    elseif x[i] < -sides[i]/2 
      x[i] = x[i] + sides[i]
    end
    x[i] = x[i] + center[i]
  end
end

# Without modifying input x

@inline function wrapone(x :: AbstractVector{Float64},
                         sides :: Vector{Float64},
                         center :: AbstractVector{Float64})
  xnew = copy(x)
  wrapone!(xnew,sides,center)
  return xnew
end
