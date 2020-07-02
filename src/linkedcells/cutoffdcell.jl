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
                      d_in_cutoff :: CutoffDistances)

  # Find the 1D index of the cell corresponding to the given 3D indexes
  jcell = icell3D_periodic(box.nc,i,j,k)

  # Check which is the index in the cell vector that is associated to the
  # input cell to be considered 
  index_cell_vector = findfirst( x -> x == jcell, lc_solvent.cell)
  if index_cell_vector == nothing
    return
  end

  # Cycle of the atoms of the solvent in this cell, computing the distances
  # and annotating the distances and the atoms of those smaller than the cutoff
  jat = lc_solvent.firstatom[index_cell_vector]
  while jat > 0
    yat = @view(x_solvent[jat,1:3])
    d = distance(box,xat,yat)
    if d <= box.cutoff
      d_in_cutoff.nd[1] += 1
      # If the number of distances found is greater than maxdim,
      # we need to increase the size of the vectors by 50%
      maxdim = length(d_in_cutoff.d)
      if d_in_cutoff.nd[1] > maxdim
        resize!(d_in_cutoff.d,round(Int64,round(Int64,1.5*maxdim)))
        resize!(d_in_cutoff.iat,round(Int64,round(Int64,1.5*maxdim)))
        resize!(d_in_cutoff.jat,round(Int64,round(Int64,1.5*maxdim)))
        resize!(d_in_cutoff.imol,round(Int64,round(Int64,1.5*maxdim)))
        resize!(d_in_cutoff.jmol,round(Int64,round(Int64,1.5*maxdim)))
        maxdim = length(d_in_cutoff.d)
        @. d_in_cutoff.d[d_in_cutoff.nd[1]+1:maxdim] = 0.
        @. d_in_cutoff.iat[d_in_cutoff.nd[1]+1:maxdim] = 0
        @. d_in_cutoff.jat[d_in_cutoff.nd[1]+1:maxdim] = 0
        @. d_in_cutoff.imol[d_in_cutoff.nd[1]+1:maxdim] = 0
        @. d_in_cutoff.jmol[d_in_cutoff.nd[1]+1:maxdim] = 0
      end
      d_in_cutoff.iat[d_in_cutoff.nd[1]] = iat
      d_in_cutoff.jat[d_in_cutoff.nd[1]] = jat
      d_in_cutoff.d[d_in_cutoff.nd[1]] = d
    end
    jat = lc_solvent.nextatom[jat]
  end

end
