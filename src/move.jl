#
# Function that generates a new random position for a molecule
#
# the new position is returned in x, a previosly allocated array
#

struct MoveAux
  oldcm :: Vector{Float64}
  newcm :: Vector{Float64}
  angles :: Vector{Float64}
  A :: Matrix{Float64}
end
MoveAux() = MoveAux(zeros(3), zeros(3), zeros(3), zeros(3,3))

function random_move!(jfmol :: Int64, jlmol :: Int64, x_solvent :: Array{Float64},
                      sizes :: Vector{Float64}, solute_center :: Vector{Float64}, 
                      x_solvent_random :: Array{Float64}, aux :: MoveAux )

  # Generate random coordiantes for the center of mass
  @. aux.newcm = -sizes/2 + rand(Float64)*sizes + solute_center 

  # Generate random rotation angles 
  @. aux.angles = (2*pi)*rand(Float64)

  # Copy the coordinates of the molecule chosen to the random-coordinates vector
  iatom = 0
  for i in jlmol:jfmol 
    iatom = iatom + 1
    x_solvent_random[iatom,1] = x_solvent[i,1]
    y_solvent_random[iatom,2] = x_solvent[i,2]
    z_solvent_random[iatom,3] = x_solvent[i,3]
  end

  # Move molecule to new position
  move!(x_solvent_random, aux)

end

#

#
# function move: Translates and rotates a molecule according
# to the desired input center of coordinates and Euler rotations
# modifyies the vector x
#

function move!(x :: Array{Float64}, aux :: MoveAux)
  
  n = length(x)
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





