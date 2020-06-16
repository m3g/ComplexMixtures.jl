#
# Function that generates a new random position for a molecule
#
# the new position is returned in x, a previosly allocated array
#

function random_move!(jfmol :: Int64, jlmol :: Int64, x_solvent :: Array{Float64},
                      sides :: Vector{Float64}, solute_center :: Vector{Float64}, 
                      x_solvent_random :: Array{Float64}, aux :: MoveAux )

  # To avoid boundary problems, the center of coordinates are generated in a 
  # much larger region, and wrapped aftwerwards
  scale = 100.

  # Generate random coordiantes for the center of mass
  @. aux.newcm = -scale*sides/2 + rand(Float64)*scale*sides + solute_center 

  # Generate random rotation angles 
  @. aux.angles = (2*pi)*rand(Float64)

  # Copy the coordinates of the molecule chosen to the random-coordinates vector
  iatom = 0
  for i in jfmol:jlmol 
    iatom = iatom + 1
    x_solvent_random[iatom,1] = x_solvent[i,1]
    x_solvent_random[iatom,2] = x_solvent[i,2]
    x_solvent_random[iatom,3] = x_solvent[i,3]
  end

  # Move molecule to new position
  move!(x_solvent_random, aux)
  
  # Wrap coordinates relative to solute center 
  wrap!(x_solvent_random,sides,solute_center)

end

