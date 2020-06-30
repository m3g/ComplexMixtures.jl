#
# function move: Translates and rotates a molecule according
# to the desired input center of coordinates and Euler rotations
# modifyies the vector x
#

function move!(x :: AbstractArray{Float64}, aux :: MoveAux)
  
  n = size(x,1)
  centerofcoordinates!(aux.oldcm,x)
  eulermat!(aux)

  # Just to simplify the code (only name assignments)
  oldcm = aux.oldcm
  newcm = aux.newcm
  A = aux.A
  @. aux.x = x
  y = aux.x

  for i in 1:n
    y[i,1] = x[i,1] - oldcm[1]
    y[i,2] = x[i,2] - oldcm[2]
    y[i,3] = x[i,3] - oldcm[3]
    x[i,1] = newcm[1] + y[i,1]*A[1,1] + y[i,2]*A[2,1] + y[i,3]*A[3,1]    
    x[i,2] = newcm[2] + y[i,1]*A[1,2] + y[i,2]*A[2,2] + y[i,3]*A[3,2]    
    x[i,3] = newcm[3] + y[i,1]*A[1,3] + y[i,2]*A[2,3] + y[i,3]*A[3,3]    
  end

end





