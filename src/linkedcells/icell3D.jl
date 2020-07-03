#
# Returns the i,j,k coordinates of the cell from the coordinates of an atom and box properties
#

function icell3D( iat :: Int64, x :: Array{Float64}, box :: Box )
  i = trunc(Int64,(x[iat,1]-box.xmin[1])/box.cutoff)+1
  j = trunc(Int64,(x[iat,2]-box.xmin[2])/box.cutoff)+1
  k = trunc(Int64,(x[iat,3]-box.xmin[3])/box.cutoff)+1
  return i, j, k
end

icell3D(x :: Array{Float64}, box :: Box) = icell3D(1,x,box)

