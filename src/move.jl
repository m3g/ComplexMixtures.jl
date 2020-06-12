#
# function move: Translates and rotates a molecule according
# to the desired input center of coordinates and Euler rotations
# modifyies the vector x
#

function move!(x :: Array{Float64}, aux :: MoveAux)
  
  n = size(x,1)
  centerofcoordinates!(aux.oldcm,x)
  eulermat!(aux)

  # Just to simplify the code (only name assignments)
  oldcm = aux.oldcm
  newcm = aux.newcm
  A = aux.A

  for i in 1:n
    x[i,1] = x[i,1] - oldcm[1]
    x[i,2] = x[i,2] - oldcm[2]
    x[i,3] = x[i,3] - oldcm[3]
    x[i,1] = newcm[1] + x[i,1]*A[1,1] + x[i,2]*A[2,1] + x[i,3]*A[3,1]    
    x[i,2] = newcm[2] + x[i,1]*A[1,2] + x[i,2]*A[2,2] + x[i,3]*A[3,2]    
    x[i,3] = newcm[3] + x[i,1]*A[1,3] + x[i,2]*A[2,3] + x[i,3]*A[3,3]    
  end

end





