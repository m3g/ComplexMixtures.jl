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


