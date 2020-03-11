# 
# Function that fills the phatnom boxes for periodic boundary conditions
# with the linked cell method
#
# On input: nboxes[3]: is the vector containing the number of cells in each
#                     direction (except phantom ones)
#           iatomfirst: Is the first atom of each cell already filled up for
#                       real cells. 
#
# On output: iatomfirst will be filled up completely in the borders with replicas
#            of the corresponding iatomfirst parameters
#

function phantomcells!(nboxes,iatomfirst) 

  # Vertices

  iatomfirst[0,0,0] = iatomfirst[nboxes[1],nboxes[2],nboxes[3]]

  iatomfirst[nboxes[1]+1,0,0] = iatomfirst[1,nboxes[2],nboxes[3]]
  iatomfirst[0,nboxes[2]+1,0] = iatomfirst[nboxes[1],1,nboxes[3]]
  iatomfirst[0,0,nboxes[3]+1] = iatomfirst[nboxes[1],nboxes[2],1]

  iatomfirst[0,nboxes(2)+1,nboxes(3)+1] = iatomfirst(nboxes(1),1,1)
  iatomfirst[nboxes(1)+1,0,nboxes(3)+1] = iatomfirst(1,nboxes(2),1)
  iatomfirst[nboxes(1)+1,nboxes(2)+1,0] = iatomfirst(1,1,nboxes(3))

  iatomfirst[nboxes[1]+1,nboxes[2]+1,nboxes[3]+1] = iatomfirst[1,1,1]

  # Axes 

  for i in 1:nboxes[1]
    iatomfirst[i,0,0] = iatomfirst[i,nboxes[2],nboxes[3]]
    iatomfirst[i,nboxes[2]+1,0] = iatomfirst[i,1,nboxes[3]]
    iatomfirst[i,0,nboxes[3]+1] = iatomfirst[i,nboxes[2],1]
    iatomfirst[i,nboxes[2]+1,nboxes[3]+1] = iatomfirst[i,1,1]
  end
  for j in 1:nboxes[2]
    iatomfirst[0,j,0] = iatomfirst[nboxes[1],j,nboxes[3]]
    iatomfirst[nboxes[1]+1,j,0] = iatomfirst[1,j,nboxes[3]]
    iatomfirst[0,j,nboxes[3]+1] = iatomfirst[nboxes[1],j,1]
    iatomfirst[nboxes[1]+1,j,nboxes[3]+1] = iatomfirst[1,j,1]
  end
  for k in 1:nboxes[3]
    iatomfirst[0,0,k] = iatomfirst[nboxes[1],nboxes[2],k]
    iatomfirst[nboxes[1]+1,0,k] = iatomfirst[1,nboxes[2],k]
    iatomfirst[0,nboxes[2]+1,k] = iatomfirst[nboxes[1],1,k]
    iatomfirst[nboxes[1]+1,nboxes[2]+1,k] = iatomfirst[1,1,k]
  end

  # Faces

  for i in 1:nboxes[2]
    for k in 1:nboxes[3]
      iatomfirst[0,j,k] = iatomfirst[nboxes[1],j,k]
      iatomfirst[nboxes[1]+1,j,k] = iatomfirst[1,j,k]
    end
  end
  for i in 1:nboxes[1]
    for k in 1:nboxes[3]
      iatomfirst[i,0,k] = iatomfirst[i,nboxes[2],k]
      iatomfirst[i,nboxes[2]+1,k] = iatomfirst[i,1,k]
    end
  end
  for i in 1:nboxes[1]
    for j in 1:nboxes[2]
      iatomfirst[i,j,0] = iatomfirst[i,j,nboxes[3]]
      iatomfirst[i,j,nboxes[3]+1] = iatomfirst[i,j,1]
    end
  end

  return
end
