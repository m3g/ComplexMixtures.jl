#
# given the 1D representation of the index of a 3D cell, returns the 3 indexes of the
# coordinates in 3 dimensions
#

function ijkcell(n :: Vector{Int64}, i1D :: Int64)
  jk = i1D%(n[2]*n[3])
  if jk == 0
    jk = n[2]*n[3]
  end
  j, k = ijcell(n[2:3],jk)
  i = round(Int64,(i1D-jk)/(n[2]*n[3]))+1
  return i, j, k
end

# Given the 1D represetation of the index of a 2D cell, returns the 2
# indexes of the 2D representation

function ijcell(n :: Vector{Int64}, i1D :: Int64)
  j = i1D%n[2]
  if j == 0
    j = n[2]
  end
  i = round(Int64,(i1D - j)/n[2])+1
  return i, j
end

