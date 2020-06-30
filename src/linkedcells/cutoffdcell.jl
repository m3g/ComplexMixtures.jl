#
# Function that computes all distance of a point "xat" to the atoms of the solvent found in
# the linked cell corresponding to indexes i, j, and k
#
# Modifies the data of d_in_cutoff
#

function cutoffdcell!(iat :: Int64, xat :: AbstractArray{Float64},
                      x_solvent :: AbstractArray{Float64},
                      lc_solvent :: LinkedCells,
                      box :: Box,
                      i :: Int64, j :: Int64, k :: Int64,
                      d_in_cutoff :: CutoffDistances,
                      nd :: Int64)

  # Find the 1D index of the cell corresponding to the given 3D indexes
  jcell = icell3D_periodic(box.nc,i,j,k)

  # Check which is the index in the cell vector that is associated to the
  # input cell to be considered 
  index_cell_vector = my_searchsortedfirst(lc_solvent.cell,jcell)
  if index_cell_vector == 0
    return
  end

  # Cycle of the atoms of the solvent in this cell, computing the distances
  # and annotating the distances and the atoms of those smaller than the cutoff
  jat = lc.firstatom[index_cell_vector]
  while jat > 0
    yat = @view(x_solvent[jat,1:3])
    d = distance(xat,yat)
    if d <= box.cutoff
      nd = nd + 1
      # If the number of distances found is greater than maxdim,
      # we need to increase the size of the vectors by 10%
      maxdim = length(d_in_cutoff.d)
      if nd > maxdim
        resize!(d_in_cutoff.d,round(Int64,round(Int64,1.1*maxdim)))
        resize!(d_in_cutoff.iat,round(Int64,round(Int64,1.1*maxdim)))
        resize!(d_in_cutoff.jat,round(Int64,round(Int64,1.1*maxdim)))
      end
      d_in_cutoff[nd] = iat
      d_in_cutoff[nd] = jat
      d_in_cutoff.d[nd] = d
    end
    jat = lc_solvent.nextatom[jat]
  end

end
