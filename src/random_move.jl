#
# Function that generates a new random position for a molecule
#
# the new position is returned in x, a previously allocated array
#
# x_solvent_random might be a view of the array that contains all the solvent
# molecules
#

function random_move!(x_ref :: AbstractArray{Float64},
                      irefatom :: Int64, sides :: Vector{Float64},
                      x_new :: AbstractArray{Float64}, aux :: MoveAux )

  # To avoid boundary problems, the center of coordinates are generated in a 
  # much larger region, and wrapped aftwerwards
  scale = 100.

  # Generate random coordiantes for the center of mass
  @. aux.newcm = -scale*sides/2 + rand(Float64)*scale*sides

  # Generate random rotation angles 
  @. aux.angles = (2*pi)*rand(Float64)

  # Copy the coordinates of the molecule chosen to the random-coordinates vector
  for i in 1:size(x_new,1)
    x_new[i,1] = x_ref[i,1]
    x_new[i,2] = x_ref[i,2]
    x_new[i,3] = x_ref[i,3]
  end
  
  # Take care that this molecule is not split by periodic boundary conditions, by
  # wrapping its coordinates around its reference atom
  wrap!(x_new,sides,@view(x_ref[irefatom,1:3]))

  # Move molecule to new position
  move!(x_new,aux)

end

