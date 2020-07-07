#
# Returns the i,j,k coordinates of the cell from the coordinates of an atom and box properties
#

function icell3D(x :: AbstractVector{Float64}, box)
  wrapone!(x,box.sides)
  i = trunc(Int64,(x[1]+box.sides[1])/box.cutoff)+1
  j = trunc(Int64,(x[2]+box.sides[2])/box.cutoff)+1
  k = trunc(Int64,(x[3]+box.sides[3])/box.cutoff)+1
  return i, j, k
end


