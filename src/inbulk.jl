#
# Function that returns if a distance is in the bulk region or not, according
# to the options
#
function inbulk(d :: Float64, R :: Result)
  if R.options.usecutoff
    if d >= R.dbulk && d < R.cutoff 
      return true
    end
  else
    if d >= R.dbulk
      return true
    end
  end
  return false
end
