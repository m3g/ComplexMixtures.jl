#
# Function that returns the first atom in the linked list of a cell,
# given the indexes of the cell in 3D
#

# From the 1D index of the cell
function first_atom_in_cell(icell :: Int64, lc :: LinkedCells) 
  index_atom = findfirst( i -> i == icell, lc.cell )
  return lc.firstatom(index_atom)
end

# From the 3D indexes of the cell
function first_atom_in_cell(i :: Int64,j :: Int64, k :: Int64, box :: Box, lc :: LinkedCells)
  icell = icell3D(box.nc,i,j,k)
  return first_atom_in_cell( icell, lc )
end

# Directly from the coordinates of a given atom
function first_atom_in_cell( x :: AbstractArray, box :: Box, lc :: LinkedCells )
  icell = icell3D(x, box)
  return first_atom_in_cell( icell, lc )
end

