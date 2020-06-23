#
# firstatom is a vector that contains for each index i3D = index3D(i,j,k) of a box, which
# is the first atom of the list of atoms in that tox
#
# nextatom is a vector that contains the index of the next atom of that box given the index
# of the previous atom (starting with firstatom).
#
# firstatom is filled up by to which box each atom pertains. The index of the box is saved
# in vector ibox, and these two vectors are ordered according to the ibox indices (from 
# greater to smaller).  
#
# Next, the loop of computations is perfomed on the ibox vector, cyling over repeated 
#

struct LinkedCells

  cell :: Vector{Int64}
  firstatom :: Vector{Int64}
  nextatom :: Vector{Int64}
  
  # Auxiliary arrays
  xmin :: Vector{Float64}
  xmax :: Vector{Float64}
  nc :: Vector{Int64}
  iaux :: Vector{Int64}
  xaux :: Vector{Int64}
  yaux :: Vector{Int64}

end

LinkedCells(n) = LinkedCells( zeros(Int64,n), # cell
                              zeros(Int64,n), # firstatom
                              zeros(Int64,n), # nextatom
                              zeros(3), # xmin
                              zeros(3), # xmax
                              zeros(Int64,3), # nc
                              zeros(Int64,n), # iaux
                              zeros(Int64,n), # xaux
                              zeros(Int64,n)  # yaux
                            )

