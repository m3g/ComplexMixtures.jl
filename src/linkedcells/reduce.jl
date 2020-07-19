#
# Reduce dc structure in the case of parallel computations
#

function reduce!(dc :: Vector{CutoffDistances})
  nthreads = Threads.nthreads()
  # Total number of cutoff distances
  ndtot = 0
  for ithread in 1:nthreads
    ndtot += dc[ithread].nd[1]
  end
  if ndtot > dc[1].maxdim[1]
    increase_size!(dc[1],round(Int64,1.5*ndtot))
  end
  for ithread in 2:nthreads
    for i in 1:dc[ithread].nd[1]
      j = dc[1].nd[1] + i
      dc[1].d[j] = dc[ithread].d[i]
      dc[1].iat[j] = dc[ithread].iat[i]
      dc[1].jat[j] = dc[ithread].jat[i]
      dc[1].imol[j] = dc[ithread].imol[i]
      dc[1].jmol[j] = dc[ithread].jmol[i]
    end
    dc[1].nd[1] += dc[ithread].nd[1]
  end
end

# If there was only a single thread, just return

function reduce!(dc :: CutoffDistances)
  return
end
























