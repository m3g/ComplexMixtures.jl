#
# firstatom is a vector that contains for each index i3D = index3D(i,j,k) of a box, which
# is the first atom of the list of atoms in that tox
#
# nextatom is a vector that contains the index of the next atom of that box given the index
# of the previous atom (starting with firstatom).
#

struct LinkedCells
  firstatom :: Vector{Int64}
  nextatom :: Vector{Int64}
end

LinkedCells(n) = LinkedCells( zeros(Int64,n), # firstatom (actual required size is nc1*nc2*nc*3
                              zeros(Int64,n), # nextatom
                            )

