#
# Function that computes all distance of a point "xat" to the atoms of the solvent found in
# the linked cell corresponding to indexes i, j, and k
#
# Modifies the data of d_in_cutoff

function cutoffdcell!(xat :: AbstractArray,
                      x_solvent :: Array,
                      lc_solvent :: LinkedCells,
                      box :: Box,
                      i :: Int64, j :: Int64, k :: Int64,
                      d_in_cutoff :: CutoffDistances)

  jcell = icell3D(box.nc,i,j,k)
  index_cell_vector = findfirst( j -> j == jcell, lc_solvent.cell ) 
  jat = lc.firstatom[index_cell_vector]
  while jat > 0
    yat = @view(x_solvent[jat,1:3])
    d = distance(xat,yat)
    

    jat = lc.nextatom[jat]
  end
    
  while igroup2 != 0 
    jj = data.groups.group2(igroup2)

    if ii == jj
      igroup2 = data.lists.iatomnext(igroup2) 
      continue
    end

    x, y, z = movephantomcoor(data,jj,ibox,jbox,kbox)

    d2 = ( data.frame.x[ii] - x )^2 + ( data.frame.y[ii] - y )^2 + ( data.frame.z[ii] - z )^2

    if d2 < data.lists.cutoff2 
      data.smalld.n = data.smalld.n + 1
      if data.smalld.n > data.smalld.nmax
        memerror = true
      else
        data.smalld.index[data.smalld.n,1] = igroup1
        data.smalld.index[data.smalld.n,2] = igroup2
        data.smalld.d[data.smalld.n] = sqrt(d2)
      end
    end

    igroup2 = data.iatomnext(igroup2)
  end

  return memerror

end
