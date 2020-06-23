#
# Given the 3D indexes of a cell, returns the index of the cell in a 
# continuous vector
#
icell3D(nc,i,j,k) = (i-1)*nc[2]*nc[3] + (j-1)*nc[3] + k

