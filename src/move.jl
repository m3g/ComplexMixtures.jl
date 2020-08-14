#
# function move: Translates and rotates a molecule according
# to the desired input center of coordinates and Euler rotations
# modifyies the vector x
#

function move!(x :: AbstractArray{Float64}, aux :: MoveAux)
  
  # Compute center of coordinates of x
  centerofcoordinates!(aux.oldcm,x)
 
  # Obtain rotation matrix
  eulermat!(aux)

  # Number of atoms
  nx = size(x,2)

  # Just to simplify the code (only name assignments)
  oldcm = aux.oldcm
  newcm = aux.newcm
  A = aux.A
  for i in 1:nx
    aux.x[1,i] = x[1,i]
    aux.x[2,i] = x[2,i]
    aux.x[3,i] = x[3,i]
  end
  y = aux.x

  for i in 1:nx
    y[1,i] = x[1,i] - oldcm[1]
    y[2,i] = x[2,i] - oldcm[2]
    y[3,i] = x[3,i] - oldcm[3]
    x[1,i] = newcm[1] + y[1,i]*A[1,1] + y[2,i]*A[1,2] + y[3,i]*A[1,3]    
    x[2,i] = newcm[2] + y[1,i]*A[2,1] + y[2,i]*A[2,2] + y[3,i]*A[2,3]    
    x[3,i] = newcm[3] + y[1,i]*A[3,1] + y[2,i]*A[3,2] + y[3,i]*A[3,3]    
  end

  return nothing
end

