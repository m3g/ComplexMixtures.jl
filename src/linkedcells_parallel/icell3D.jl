#
# Returns the i,j,k coordinates of the cell from the coordinates of an atom and box properties
#

function icell3D(x :: AbstractVector{Float64}, box)
  i = trunc(Int64,(x[1]+box.sides[1]/2)/box.l[1])+1
  j = trunc(Int64,(x[2]+box.sides[2]/2)/box.l[2])+1
  k = trunc(Int64,(x[3]+box.sides[3]/2)/box.l[3])+1
  return i, j, k
end

#
# given the 1D representation of the index of a 3D cell, returns the 3 indexes of the
# coordinates in 3 dimensions
#

# if the number of cells in each dimension is different

function icell3D(n :: Vector{Int64}, i1D :: Int64)
  jk = i1D%(n[2]*n[3])
  if jk == 0
    jk = n[2]*n[3]
  end
  j, k = icell2D(n[2:3],jk)
  i = round(Int64,(i1D-jk)/(n[2]*n[3]))+1
  return i, j, k
end

# If the number of cells in each dimension is the same

function icell3D(n :: Int64, i1D :: Int64)
  n2 = n^2
  jk = i1D%n2
  if jk == 0
    jk = n2
  end
  j, k = icell2D(n,jk)
  i = round(Int64,(i1D-jk)/(n2))+1
  return i, j, k
end

# Given the 1D represetation of the index of a 2D cell, returns the 2
# indexes of the 2D representation

# If the two dimensions are different

function icell2D(n :: Vector{Int64}, i1D :: Int64)
  j = i1D%n[2]
  if j == 0
    j = n[2]
  end
  i = round(Int64,(i1D - j)/n[2])+1
  return i, j
end

# If the two dimensions are the same

function icell2D(n :: Int64, i1D :: Int64)
  j = i1D%n
  if j == 0
    j = n
  end
  i = round(Int64,(i1D - j)/n)+1
  return i, j
end



