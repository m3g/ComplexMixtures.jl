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

# Special case in which the cells are periodic: if the indexes received are outside the
# boundaries (being 0 or n+1), return the indexes of the corresponding periodic cell on the
# other side

function icell3D_periodic(nc, i, j, k)

  # Vertices
  i == 0         && j == 0         && k == 0         && return icell3D(nc, nc[1], nc[2], nc[3])

  i == nc[1] + 1 && j == 0         && k == 0         && return icell3D(nc,   1  , nc[2], nc[3])
  i == 0         && j == nc[2] + 1 && k == 0         && return icell3D(nc, nc[1],   1  , nc[3])
  i == 0         && j == 0         && k == nc[3] + 1 && return icell3D(nc, nc[1], nc[2],   1  )

  i == 0         && j == nc[2] + 1 && k == nc[3] + 1 && return icell3D(nc, nc[1],   1  ,   1  )
  i == nc[1] + 1 && j == 0         && k == nc[3] + 1 && return icell3D(nc,   1  , nc[2],   1  )
  i == nc[1] + 1 && j == nc[2] + 1 && k == 0         && return icell3D(nc,   1      1  , nc[3])

  i == nc[1] + 1 && j == nc[2] + 1 && k == nc[3] + 1 && return icell3D(nc,   1  ,   1,     1  )

  # Axes 
  j == 0         && k == 0         && return icell3D(nc,   i  , nc[2], nc[3])
  j == nc[2] + 1 && k == 0         && return icell3D(nc,   i  ,   1  , nc[3])
  j == 0         && k == nc[3] + 1 && return icell3D(nc,   i  , nc[2],   1  )
  j == nc[2] + 1 && k == nc[3] + 1 && return icell3D(nc,   i  ,   1  ,   1  )

  i == 0         && k == 0         && return icell3D(nc, nc[1],  j  , nc[3])
  i == nc[1] + 1 && k == 0         && return icell3D(nc,   1  ,  j  , nc[3])
  i == 0         && k == nc[3] + 1 && return icell3D(nc, nc[1],  j  ,   1  )
  i == nc[1] + 1 && k == nc[3] + 1 && return icell3D(nc,   1  ,  j  ,   1  )

  i == 0         && j == 0         && return icell3D(nc, nc[1], nc[2],  k  )
  i == nc[1] + 1 && j == 0         && return icell3D(nc,   1  , nc[2],  k  )
  i == 0         && j == nc[2] + 1 && return icell3D(nc, nc[1],   1  ,  k  )
  i == nc[1] + 1 && j == nc[2] + 1 && return icell3D(nc,   1  ,   1  ,  k  )

  # Faces
  i == 0         && return icell3D(nc, nc[1],  j   ,   k  )
  i == nc[1] + 1 && return icell3D(nc,   1  ,  j   ,   k  )

  j == 0         && return icell3D(nc,   i  , nc[2],   k  )
  j == nc[2] + 1 && return icell3D(nc,   i  ,   1  ,   k  )

  k == 0         && return icell3D(nc,   i  ,   j  , nc[3])
  k == nc[3] + 1 && return icell3D(nc,   i  ,   j  ,   1  )

  # not outside, return the indexes of a regular cell
  return icell3D(nc,i,j,k)

end
