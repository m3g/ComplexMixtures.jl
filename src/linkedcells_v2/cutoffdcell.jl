#
# Function that computes all distance of a point "xat" to the atoms of the solvent found in
# the linked cell corresponding to indexes i, j, and k
#
# Modifies the data of dc
#

function cutoffdcell!(iat :: Int64, xat :: AbstractArray{Float64},
                      x_solvent :: AbstractArray{Float64},
                      lc_solvent :: LinkedCells,
                      box :: Box,
                      i :: Int64, j :: Int64, k :: Int64,
                      dc :: CutoffDistances)

  # Check if this box needs to be wrapped. If so, the distance calculation has to take
  # that in consideration
  wrapped = false
  if i == 0 || j == 0 || k == 0 || i == box.nc[1]+1 || j == box.nc[2]+1 || k == box.nc[3]+1
    wrapped = true
  end

  # Get the indexes of the current cell considering the possible wrap associated
  # to periodic boundary conditions
  i, j, k = wrap_cell(box.nc,i,j,k)
  icell = icell1D(box.nc,i,j,k)

  # Cycle of the atoms of the solvent in this cell, computing the distances
  # and annotating the distances and the atoms of those smaller than the cutoff
  jat = lc_solvent.firstatom[icell]
  while jat > 0
    if wrapped 
      d = distance(box.sides,@view(x_solvent[jat,1:3]),xat)
    else
      d = distance(@view(x_solvent[jat,1:3]),xat)
    end
    if d <= box.cutoff
      dc.nd[1] += 1
      # If the number of distances found is greater than maxdim,
      # we need to increase the size of the vectors by 50%
      if dc.nd[1] > dc.maxdim[1]
        increase_size!(dc)
      end
      dc.iat[dc.nd[1]] = iat
      dc.jat[dc.nd[1]] = jat
      dc.d[dc.nd[1]] = d
    end
    jat = lc_solvent.nextatom[jat]
  end

end
