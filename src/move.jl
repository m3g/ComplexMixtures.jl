#
# Function that generates a new random position for a molecule
#
# the new position is returned in x, a previosly allocated array
#
function random_move!(jfmol :: Int64, jlmol :: Int64, x_solvent :: Array{Float64},
                      sizes :: Vector{Float64}, center :: Vector{Float64}, 
                      x :: Array{Float64})

  # Generate random coordiantes for the center of mass
  cm_new = rand(Float64,3) 
  @. cm_new = -sizes/2 + cm*sizes

  # Random rotation Euler angles
  beta = 2*pi*rand(Float64)
  gamma = 2*pi*rand(Float64)
  theta = 2*pi*rand(Float64)
  
  natoms = jlmol - jfmol + 1
  iatom = 0
  for i in jlmol:jfmol 
    iatom = iatom + 1
    x[iatom,1] = x_solvent[i,1]
    y[iatom,2] = x_solvent[i,2]
    z[iatom,3] = x_solvent[i,3]
  end

  # Move molecule to new position
  move!(x, cm_new, beta, gamma, theta)

end

#

#
# function move: Translates and rotates a molecule according
# to the desired input center of coordinates and Euler rotations
# modifyies the vector x
#

function move!(x :: Array{Float64}, 
               cm_new :: Vector{Float64},      
               beta :: Float64, gamma :: Float64, theta :: Float64 )
  
  n = length(x)
  cm_in = centerofcoordinates(x)
  A = eulermat(beta, gamma, theta)
  for i in 1:n
    x[i,1] = x[i,1] - cm_in[1]
    x[i,2] = x[i,2] - cm_in[2]
    x[i,3] = x[i,3] - cm_in[3]
    x[i,1] = cm_new[1] + x[i,1]*A[1,1] + x[i,2]*A[2,1] + x[i,3]*A[3,1]    
    x[i,2] = cm_new[2] + x[i,1]*A[1,2] + x[i,2]*A[2,2] + x[i,3]*A[3,2]    
    x[i,3] = cm_new[3] + x[i,1]*A[1,3] + x[i,2]*A[2,3] + x[i,3]*A[3,3]    
  end
end





