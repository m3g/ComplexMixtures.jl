#
# Function that initializes the linked cells by computing to each cell each atom
# belongs and filling up the firstatom and nexatom arrays.
#
# Modifies the data in the lc structure
#
                  
function initcells!(x :: Array{Float64}, box :: Box, lc :: LinkedCells )
  ifmol = 1
  ilmol = size(x,1)
  return initcells!(ifmol,ilmol,x,box,lc)
end

function initcells!(ifmol :: Int64, ilmol :: Int64, x :: Array{Float64}, box :: Box, lc :: LinkedCells )

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
  for iat in ifmol:ilmol
    ic, jc, kc = icell3D(iat,x,box)
    icell = icell1D(box.nc,ic,jc,kc)
    lc.nextatom[iat] = lc.firstatom[icell]
    lc.firstatom[icell] = iat
  end

end

 
