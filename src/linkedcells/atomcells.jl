#
# Given the vector of coordinates and the cutoff (which corresponds to the
# cell side) returns the index of the cell to which each atom belongs
#
# We are assuming here that the coordinates are already wrapped around the center
# of coordinates of interest, in such a way that periodic boundary conditions
# are not a problem here
#
function atomcells(n :: Int64, x :: Array{Float64}, cutoff :: Float64, 
                   cells :: Vector{Int64}) 

  cmin = Vector{Float64}(undef,3)
  for i in 1:3
    cmin[i] = minimum(x[:,i])
  end

  for i in 1:n
    ix = trunc((x[i,1]-cmin[1])/cutoff)+1
    iy = trunc((x[i,2]-cmin[2])/cutoff)+1
    iz = trunc((x[i,3]-cmin[3])/cutoff)+1
    cells[i] = icell3D(ix,iy,iz)
  end

end
