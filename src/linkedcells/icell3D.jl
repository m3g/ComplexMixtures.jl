#
# returns the index of the linked cell, in the 1D representation, to which an atom belongs
#

# From the indexes of the cell
icell3D(nc :: Vector{Int64}, i :: Int64, j :: Int64, k :: Int64) = (i-1)*nc[2]*nc[3] + (j-1)*nc[3] + k

# From the coordinates of an atom
function icell3D( x :: AbstractArray, cutoff :: Float64, lc :: LinkedCells )
  i = trunc(Int64,(x[1]-lc.xmin[1])/cutoff)+1
  j = trunc(Int64,(x[2]-lc.xmin[2])/cutoff)+1
  k = trunc(Int64,(x[3]-lc.xmin[3])/cutoff)+1
  return icell3D(lc.nc,i,j,k)
end

