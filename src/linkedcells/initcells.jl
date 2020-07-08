#
# Function that initializes the linked cells by computing to each cell each atom
# belongs and filling up the firstatom and nexatom arrays.
#
# Modifies the data in the lc structure
#
                  
function initcells!(x :: AbstractArray{Float64}, box :: Box, lc :: LinkedCells )

  # Wrap coordinates relative to the origin. This is necessary to speedup calculations,
  # as the wrapping of coordinates to compute distanes will be only necessary for cells
  # at the borders of the simulation box. Not having the compute minimum images between
  # pairs of atoms actually speeds up the calculations significantly, because the calculation
  # of pairwise distances is the most repeated calculation 
  wrap!(x,box.sides)

  nboxes = box.nc[1]*box.nc[2]*box.nc[3] 
  if length(lc.firstatom) < nboxes
    resize!(lc.firstatom,nboxes)
  end

  # Reset arrays
  for i in 1:nboxes
    lc.firstatom[i] = 0
  end
  @. lc.nextatom = 0

  # Initialize cell, firstatom and nexatom
  for iat in 1:size(x,1)
    ic, jc, kc = icell3D(@view(x[iat,1:3]),box)
    icell = icell1D(box.nc,ic,jc,kc)
    lc.nextatom[iat] = lc.firstatom[icell]
    lc.firstatom[icell] = iat
  end

end

 
