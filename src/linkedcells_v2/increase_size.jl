#
# Increases size of arrays in dc structure if needed
#
function increase_size!(dc :: CutoffDistances)
  dc.maxdim[1] = round(Int64,1.5*dc.maxdim[1])
  resize!(dc.d,dc.maxdim[1])
  resize!(dc.iat,dc.maxdim[1])
  resize!(dc.jat,dc.maxdim[1])
  resize!(dc.imol,dc.maxdim[1])
  resize!(dc.jmol,dc.maxdim[1])
  @. dc.d[dc.nd[1]+1:dc.maxdim[1]] = 0.
  @. dc.iat[dc.nd[1]+1:dc.maxdim[1]] = 0
  @. dc.jat[dc.nd[1]+1:dc.maxdim[1]] = 0
  @. dc.imol[dc.nd[1]+1:dc.maxdim[1]] = 0
  @. dc.jmol[dc.nd[1]+1:dc.maxdim[1]] = 0
end

