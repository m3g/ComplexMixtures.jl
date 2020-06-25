#
# returns the index of the linked cell, in the 1D representation, to which an atom belongs
#

# From the indexes of the cell
icell3D(nc :: Vector{Int64}, i :: Int64, j :: Int64, k :: Int64) = (i-1)*nc[2]*nc[3] + (j-1)*nc[3] + k

# From the coordinates of an atom and box properties
function icell3D( x :: AbstractArray, box :: Box )
  i = trunc(Int64,(x[1]-box.xmin[1])/box.cutoff)+1
  j = trunc(Int64,(x[2]-box.xmin[2])/box.cutoff)+1
  k = trunc(Int64,(x[3]-box.xmin[3])/box.cutoff)+1
  return icell3D(box.nc,i,j,k)
end

