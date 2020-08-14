#
# Functions that wrap the coordinates (x{N,3} array) to obtain minimum images
# around a defined center
# 
# It modifies the coordinates input vector
#


function wrap!(x :: AbstractArray{Float64}, 
               sides :: AbstractVector{Float64}, 
               center :: AbstractVector{Float64})
  for i in 1:size(x,2)
    wrapone!(@view(x[1:3,i]),sides,center)
  end
  return nothing
end

@inline function wrapone!(x :: AbstractVector{Float64}, 
                          sides :: AbstractVector{Float64}, 
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
  return nothing
end

# Without modifying input x

@inline function wrapone(x :: AbstractVector{Float64},
                         sides :: AbstractVector{Float64},
                         center :: AbstractVector{Float64})
  xnew = copy(x)
  wrapone!(xnew,sides,center)
  return xnew
end

# If only the sides are provided, wrap to origin

function wrap!(x :: AbstractArray{Float64}, 
               sides :: AbstractVector{Float64})
  for i in 1:size(x,2)
    wrapone!(@view(x[1:3,i]),sides)
  end
  return nothing
end

@inline function wrapone!(x :: AbstractVector{Float64}, 
                          sides :: AbstractVector{Float64})
  for i in 1:3
    x[i] = x[i]%sides[i]
    if x[i] > sides[i]/2  
      x[i] = x[i] - sides[i] 
    elseif x[i] < -sides[i]/2 
      x[i] = x[i] + sides[i]
    end
  end
  return nothing
end

