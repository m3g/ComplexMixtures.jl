#
# The cutoffdistances function returns every distance that is smaller than the cutoff. However,
# we only want the minimum distance between the solute molecule and the solvent atoms. Therefore,
# if two distances correspond to different solute atoms but to the same solevent atom, the greater
# of those distances must be eliminiated. 
#
# This is a first step towards cleaning the list to keep only the actual minimum distance between
# any solute and solvent atoms, which is used to keep the reference atom distances as well in 
# the list for the computatio of the RDF and of the monte-carlo integration 
#

# This function will return the d_in_cutoff data modified such that only one distance corresponding
# to a solvent atom and a solute molecule is kept in the list. Modifies everything in the CutoffDistances
# structure, including nd 

function keepunique!(d_in_cutoff :: CutoffDistances, solute)

  # Just use simpler names
  nd = d_in_cutoff.nd
  iat = d_in_cutoff.iat
  jat = d_in_cutoff.jat
  imol = d_in_cutoff.imol
  jmol = d_in_cutoff.jmol
  d = d_in_cutoff.d

  # Check to which molecule of the solute each atom of the list belongs
  for i in 1:nd[1]
    imol[i] = solute.imol[iat[i]]
  end

  for i in 1:nd[1]-1
    # If this position is already a repeated one, lets not waste time
    if imol[i] == 0
      continue
    end
    for j in i+1:nd[1]
      # If this position is already a repeated one, lets not waste time
      if imol[j] == 0
        continue  
      end
      # If the solute molecules involved are the same, keep the smallest distance
      # in the smaller index, and eliminates the other
      if imol[j] == imol[i] && jat[j] == jat[i] 
        if d[j] < d[i]
          iat[i] = iat[j]
          jat[i] = jat[j]
          d[i] = d[j]
        end
        iat[j] = 0
        jat[j] = 0
        imol[j] = 0
        jmol[j] = 0
        d[j] = 0.
      end
    end
  end

  # Put the zeros at the end
  nonzero = 0
  for i in 1:nd[1]
    if imol[i] > 0 
      nonzero = nonzero + 1
      if nonzero < i
        d[nonzero] = d[i]
        iat[nonzero] = iat[i]
        jat[nonzero] = jat[i]
        imol[nonzero] = imol[i]
        jmol[nonzero] = jmol[i]
        d[i] = 0.
        iat[i] = 0
        jat[i] = 0
        imol[i] = 0
        jmol[i] = 0
      end
    end
  end

  # Update the number of distances in the list
  nd[1] = nonzero

end

