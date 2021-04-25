#
# Functions that wrap the coordinates 
# 
# It modifies the coordinates of the input vector
#

#
# Wrap to a given center of coordinates
#
function wrap!(x::AbstractVector{T}, 
               sides::T, 
               center::T) where T <: AbstractVector
  for i in eachindex(x)
    x[i] = wrapone(x[i],sides,center)
  end
  nothing
end

@inline function wrapone(x::T, sides::T, center::T) where T <: AbstractVector
  s = @. (x-center)%sides
  s = @. wrapx(s,sides) + center
  return s
end

@inline function wrapx(x,s)
  if x > s/2
    x = x - s
  elseif x < -s/2
    x = x + s
  end
  x
end

#
# Wrap to origin
#
function wrap!(x::AbstractVector{T}, sides::T) where T <: AbstractVector
  for i in eachindex(x)
    x[i] = wrapone(x[i],sides)
  end
  nothing
end

@inline function wrapone(x::T, sides::T) where T <: AbstractVector
  s = @. x%sides
  s = @. wrapx(s,sides)
  return s
end
