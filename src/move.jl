#
# function move: Translates and rotates a molecule according
# to the desired input center of coordinates and Euler rotations
# modifyies the vector x
#

function move!(x :: Array{Float64}, aux :: MoveAux)
  first = 1
  last = size(x,1)
  move!(first,last,x,aux)
end

function move!(first :: Int64, last :: Int64, x :: Array{Float64}, aux :: MoveAux)
  
  centerofcoordinates!(aux.oldcm,first,last,x)
  eulermat!(aux)

  # Just to simplify the code (only name assignments)
  oldcm = aux.oldcm
  newcm = aux.newcm
  A = aux.A
  i = 0
  for iat in first:last
    i = i + 1
    aux.x[i,1] = x[iat,1]
    aux.x[i,2] = x[iat,2]
    aux.x[i,3] = x[iat,3]
  end
  y = aux.x

  i = 0
  for iat in first:last
    i = i + 1
    y[i,1] = x[i,1] - oldcm[1]
    y[i,2] = x[i,2] - oldcm[2]
    y[i,3] = x[i,3] - oldcm[3]
    x[iat,1] = newcm[1] + y[i,1]*A[1,1] + y[i,2]*A[2,1] + y[i,3]*A[3,1]    
    x[iat,2] = newcm[2] + y[i,1]*A[1,2] + y[i,2]*A[2,2] + y[i,3]*A[3,2]    
    x[iat,3] = newcm[3] + y[i,1]*A[1,3] + y[i,2]*A[2,3] + y[i,3]*A[3,3]    
  end

end





