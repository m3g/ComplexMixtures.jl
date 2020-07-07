# Given the indexes of a cell, return the periodic cell which correspondst to
# it, if the cell is outside the main box

function wrap_cell(nc :: Vector{Int64}, i :: Int64, j :: Int64, k :: Int64)

  if i == 0
    i = nc[1]
  elseif i == nc[1] + 1
    i = 1
  end
 
  if j == 0
    j = nc[2]
  elseif j == nc[2] + 1
    j = 1
  end

  if k == 0
    k = nc[3]
  elseif k == nc[3] + 1
    k = 1
  end

  return i, j, k

end

