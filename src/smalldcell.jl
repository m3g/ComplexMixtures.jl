#
# Function that checks the distance criterium and build the list of distances
# smaller than the cutoff
#

# Modifies the data of "data.smalld" structure

function smalldcell!(data :: DistanceData, ii, igroup1, ibox, jbox, kbox)
    
  data.smalld.memerror = false

  igroup2 = data.lists.iatomfirst(ibox,jbox,kbox)
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
        data.smalld.memerror = true
      else
        data.smalld.index[data.smalld.n,1] = igroup1
        data.smalld.index[data.smalld.n,2] = igroup2
        data.smalld.d[data.smalld.n] = sqrt(d2)
      end
    end

    igroup2 = data.iatomnext(igroup2)
  end

end
