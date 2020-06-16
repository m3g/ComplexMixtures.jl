
# Computes the center of coordinates of a vector coor, and returns it in the vector
# cm, which must be provided previously allocated


# With explicit passing of indexes of a subarray

function centerofcoordinates!(cm :: Vector{Float64}, ifirst, ilast, coor :: Array{Float64})
  @. cm = 0.
  for i in ifirst:ilast
    cm[1] = cm[1] + coor[i,1]
    cm[2] = cm[2] + coor[i,2]
    cm[3] = cm[3] + coor[i,3]
  end
  @. cm = cm / (ilast-ifirst+1)
end

function centerofcoordinates!(cm :: Vector{Float64}, coor :: Array{Float64}) 
  centerofcoordinates!(cm, 1, size(coor,1), coor)
end

# without previous allocation

function centerofcoordinates(ifirst, ilast, coor :: Array{Float64} )
  cm = zeros(3)
  centerofcoordinates!(cm,ifirst,ilast,coor)
  return cm
end

function centerofcoordinates(cm :: Vector{Float64}, coor :: Array{Float64}) 
  centerofcoordinates(cm, 1, size(coor,1), coor)
end


