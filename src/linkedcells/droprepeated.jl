#
# These functions are used to list the cells containing atoms and the first
# atom of the list of atoms of each cell, without using arrays with Ncells^3
# elements.
#
# To do so, the cellatom vector, of length Natoms, is filled up with the cell
# of the atoms, with indexes listed in the firstatom vector. 
# (while the nexataom vector is filled with the following atoms
# as usual). 
#
# Then, we elimiate from the cellatom and firstatom vectors all elements
# which are repeated elements in the cellatom array, keeping only the last
# appearance of each element, which will correspond to the last "firstatom"
# found in the construction of the linked firstatom and nexatom lists. 
#
# The functions here perform this elimination of elements of the two vectors, returning
# the vectors filled up with values different from zero up to the position of the last
# unique element.
#
# On input:
# x : cell to which each atom pertains
# y : index of the atom corresponding to position in x
#
# On output:
# y : First atom of each cell, without repeated cells, put at the begining of the vector,
#     and with zeros afterwards
#

droprepeated!( lc :: LinkedCells ) =droprepeated!( lc.cell, lc.firstatom, lc.iaux, lc.xaux, lc.yaux )

function droprepeated(xin :: Vector{Int64}, yin :: Vector{Int64})
  x = copy(xin)
  y = copy(yin)
  droprepeated!(x,y)
  return y
end

function droprepeated!(x :: Vector{Int64}, y :: Vector{Int64})  
  iaux = similar(x)
  xaux = similar(x)
  yaux = similar(x)
  return droprepeated!(x,y,iaux,xaux,yaux)
end

function droprepeated!(x :: Vector{Int64}, y :: Vector{Int64}, 
                       iaux :: Vector{Int64}, xaux :: Vector{Int64}, yaux :: Vector{Int64})

  # Sort elements of x and save indexes in iaux

  n = length(x)
  for i in 1:n
    iaux[i] = i
  end
  sort!(iaux, by = iaux -> x[iaux])
  @. xaux = x[iaux]
  @. yaux = y[iaux]
  @. y = yaux

  # Remove repeated elements of x and corresponding elements in y

  i = 1
  iunique = 1
  xunique = xaux[i]
  while i < n
    i = i + 1
    while i <= n && xaux[i] == xunique
      y[iunique] = y[i]
      i = i + 1
    end
    if i < n
      iunique = iunique + 1
      xunique = xaux[i]
      y[iunique] = y[i]
    end
  end
  for i in iunique+1:n
    y[i] = 0
  end

end


