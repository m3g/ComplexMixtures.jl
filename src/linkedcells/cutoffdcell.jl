#
# Function that computes all distance of a point "x_solute" to the atoms of the solvent found in
# the linked cell corresponding to indexes i, j, and k
#
# Modifies the data of dc
#

function cutoffdcell!(iat :: Int64, x_solute :: Array{Float64},
                      x_solvent :: Array{Float64},
                      lc_solvent :: LinkedCells,
                      box :: Box,
                      i :: Int64, j :: Int64, k :: Int64,
                      dc :: CutoffDistances)

  # Get the indexes of the current cell considering the possible wrap associated
  # to periodic boundary conditions
  i, j, k, wrapped = wrap_cell(box.nc,i,j,k)
  icell = icell1D(box.nc,i,j,k)

  # Maximum number of distances stored, might be updated in the loop if required
  maxdim = length(dc.d)

  # Cycle of the atoms of the solvent in this cell, computing the distances
  # and annotating the distances and the atoms of those smaller than the cutoff
  jat = lc_solvent.firstatom[icell]
  while jat > 0
    if wrapped
      d = distance(box,x_solute,x_solvent,iat,jat)
    else
      d = distance(x_solute,x_solvent,iat,jat)
    end
    if d <= box.cutoff
      dc.nd[1] += 1
      # If the number of distances found is greater than maxdim,
      # we need to increase the size of the vectors by 50%
      if dc.nd[1] > maxdim
        maxdim = round(Int64,1.5*maxdim)
        resize!(dc.d,maxdim)
        resize!(dc.iat,maxdim)
        resize!(dc.jat,maxdim)
        resize!(dc.imol,maxdim)
        resize!(dc.jmol,maxdim)
        @. dc.d[dc.nd[1]+1:maxdim] = 0.
        @. dc.iat[dc.nd[1]+1:maxdim] = 0
        @. dc.jat[dc.nd[1]+1:maxdim] = 0
        @. dc.imol[dc.nd[1]+1:maxdim] = 0
        @. dc.jmol[dc.nd[1]+1:maxdim] = 0
      end
      dc.iat[dc.nd[1]] = iat
      dc.jat[dc.nd[1]] = jat
      dc.d[dc.nd[1]] = d
    end
    jat = lc_solvent.nextatom[jat]
  end

end
