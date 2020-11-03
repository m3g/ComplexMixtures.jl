#
# Function that reorders x vector by putting in the first positions the
# elements with values smaller than cutoff
#
function partialsort_cutoff!(x,cutoff; by = x -> x)
  iswap = 1
  for i in 1:length(x)
    if by(x[i]) < cutoff
      if iswap == i
        iswap = iswap + 1
      else  
        xtmp = x[iswap]
        x[iswap] = x[i]
        x[i] = xtmp
        iswap = iswap + 1
      end
    end
  end
  nothing
end
