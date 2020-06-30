#
# Function that initializes the linked cells by computing to each cell each atom
# belongs and filling up the firstatom and nexatom arrays.
#
# Modifies the data in the lc structure
#
                  
function initcells!(x :: AbstractArray{Float64}, box :: Box, lc :: LinkedCells )

  n = size(x,1)

  # Reset arrays
  @. lc.cell = 0
  @. lc.firstatom = 0
  @. lc.nextatom = 0

  # Compute to which cell each atom belongs
  for i in 1:n
    icell = icell3D(@view(x[i,1:3]),box)
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

 
