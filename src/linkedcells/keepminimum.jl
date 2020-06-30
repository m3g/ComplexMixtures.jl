#
# Function that resturns the structure containing the list of the distances found keeping
# only the minimum distance between pairs of molecules, the remaining of the vectors
# will be zeroed
#

# This function will return the d_in_cutoff data modified such that only the pair of atoms
# corresponding to the minimum distance between each pair of molecules is retained in the 
# arrays. All other elements will be zeroed.   

function keepminimum!(nd :: Int64, d :: Vector{Float64}, iat :: Vector{Int64}, jat :: Vector{Int64},
                      imol :: Vector{Int64}, jmol :: Vector{Int64})
  istore = 0
  for i in 1:nd-1
    if imol[iat[i]] > 0 
      istore = istore + 1 
    else
      continue
    end
    for j in i+1:nd
      if imol[iat[j]] == imol[iat[i]] && jmol[jat[j]] == jmol[jat[i]] 
        imol[iat[j]] = 0
        if d[j] < d[i]
          iat[istore] = iat[i]
          jat[istore] = jat[i]
          d[istore] = d[j]
        end
      end
    end
  end
  @. d[istore+1:nd] = 0.
  @. jat[istore+1:nd] = 0.
end

# Calling the function using the structures instead of the vectors

keepminimum!(nd,d_in_cutoff,solute,solvent) = keepminimum!(nd,d_in_cutoff.d,d_in_cutoff.iat,d_in_cutoff.jat,
                                                           solute.imol, solvent.imol)
# This is the same function, but it keeps the minimum distance between pairs of atoms, 
# independently of their molecules. Won't be used here.

function keepminimum!(nd :: Int64, d :: Vector{Float64}, iat :: Vector{Int64}, jat :: Vector{Int64})
  istore = 0
  for i in 1:nd-1
    if iat[i] > 0 
      istore = istore + 1 
    else
      continue
    end
    for j in i+1:nd
      if iat[j] == iat[i] && jat[j] == jat[i] 
        iat[j] = 0
        if d[j] < d[i]
          iat[istore] = iat[i]
          jat[istore] = jat[i]
          d[istore] = d[j]
        end
      end
    end
  end
  @. d[istore+1:nd] = 0.
  @. jat[istore+1:nd] = 0.
end

