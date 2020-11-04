#
# Returns a view of a coordinate vector corresponding to the atoms
# of a molecule with index i. n is the number of atoms of the molecule
#
function viewmol(i :: Int64, x :: Vector{T}, n :: Int64) where T <: Vf3
  first = (i-1)* n + 1
  last = first + n - 1
  return @view(x[first:last])
end

# From the structure of the solute or solvent
viewmol(i :: Int64, x :: Vector{T}, s :: Selection) where T <: Vf3 =
  viewmol(i, x, s.natomspermol)

