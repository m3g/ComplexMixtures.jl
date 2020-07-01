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
  n_cells_with_atoms = 0
  for i in 1:n
    icell = icell3D(@view(x[i,1:3]),box)
    index_icell = findfirst( ic -> ic == icell, lc.cell )
    if index_icell == nothing
      n_cells_with_atoms += 1
      index_icell = n_cells_with_atoms
      lc.cell[index_icell] = icell
    end
    lc.nextatom[i] = lc.firstatom[index_icell]
    lc.firstatom[index_icell] = i
  end

  # Remove repeated cells from cell list and first atom list
  #droprepeated!(lc)

end

 
