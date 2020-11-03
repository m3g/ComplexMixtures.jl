#
# function move: Translates and rotates a molecule according
# to the desired input center of coordinates and Euler rotations
# modifyies the vector x
#

function move!(x :: AbstractVector{T}, newcm :: T, 
               beta :: Float64, gamma :: Float64, theta :: Float64) where T <: Vf3
  
  # Compute center of coordinates of x
  cm = centerofcoordinates(x)
 
  # Obtain rotation matrix
  A = eulermat(beta, gamma, theta)

  for i in eachindex(x)
    x[i] = A*(x[i]-cm) + newcm
  end

  nothing
end




