#
# Returns a view of a coordinate vector corresponding to the atoms
# of a molecule with index i. n is the number of atoms of the molecule
#
function viewmol(i :: Int, x :: Vector{T}, n :: Int) where T
  first = (i-1)* n + 1
  last = first + n - 1
  return @view(x[first:last])
end

# From the structure of the solute or solvent
viewmol(i :: Int, x :: Vector{T}, s :: Selection) where T =
  viewmol(i, x, s.natomspermol)

