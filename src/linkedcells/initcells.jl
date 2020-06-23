#
# Function that initializes the linked cells by computing to each cell each atom
# belongs and filling up the firstatom and nexatom arrays.
#

function initcells(n, x :: Array{Float64}, cutoff, lc :: LinkedCells )

  # Reset arrays
  @. lc.cell = 0
  @. lc.firstatom = 0
  @. lc.nextatom = 0

  # Compute minimum coordinates
  for i in 1:3
    lc.xmin[i] = minimum(x[:,i])
    lc.xmax[i] = maximum(x[:,i])
    lc.nc[i] = trunc(Int64,(lc.xmax[i]-lc.xmin[i])/cutoff) + 1
  end

  # Compute to which cell each atom belongs
  for i in 1:n
    ix = trunc(Int64,(x[i,1]-lc.xmin[1])/cutoff)+1
    iy = trunc(Int64,(x[i,2]-lc.xmin[2])/cutoff)+1
    iz = trunc(Int64,(x[i,3]-lc.xmin[3])/cutoff)+1
    icell = icell3D(lc.nc,ix,iy,iz) 
    ifirst = findlast( ic -> ic == icell, lc.cell )
    if ifirst == nothing
      ifirst = i
    end
    lc.cell[i] = icell
    lc.firstatom[i] = i
    lc.nextatom[ifirst] = i
  end

  # Remove repeated cells from cell list and first atom list
  droprepeated!(lc)

end

 
