# Given the indexes of a cell, return the periodic cell which correspondst to
# it, if the cell is outside the main box

function wrap_cell(nc :: AbstractVector, i :: Int, j :: Int, k :: Int)

  if i < 1
    i = nc[1] + i
  elseif i > nc[1]
    i = i - nc[1]
  end
    
  if j < 1
    j = nc[2] + j
  elseif j > nc[2]
    j = j - nc[2]
  end
    
  if k < 1
    k = nc[3] + k
  elseif k > nc[3]
    k = k - nc[3]
  end
    
  return i, j, k

end

