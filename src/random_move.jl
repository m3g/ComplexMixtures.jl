#
# Function that generates a new random position for a molecule
#
# the new position is returned in x, a previously allocated array
#
# x_solvent_random might be a view of the array that contains all the solvent
# molecules
#

function random_move!(x_ref :: AbstractVector{T}, 
                      irefatom :: Int64,
                      sides :: T,
                      x_new :: AbstractVector{T}, RNG) where T <: Vf3

  # To avoid boundary problems, the center of coordinates are generated in a 
  # much larger region, and wrapped aftwerwards
  scale = 100.

  # Generate random coordiantes for the center of mass
  ncm1 = -scale*sides[1]/2 + random(RNG,Float64)*scale*sides[1]
  ncm2 = -scale*sides[2]/2 + random(RNG,Float64)*scale*sides[2]
  ncm3 = -scale*sides[3]/2 + random(RNG,Float64)*scale*sides[3]
  newcm = T(ncm1,ncm2,ncm3)

  # Generate random rotation angles 
  beta = (2*pi)*random(RNG,Float64)
  gamma = (2*pi)*random(RNG,Float64)
  theta = (2*pi)*random(RNG,Float64)

  # Copy the coordinates of the molecule chosen to the random-coordinates vector
  @. x_new = x_ref

  # Take care that this molecule is not split by periodic boundary conditions, by
  # wrapping its coordinates around its reference atom
  wrap!(x_new,sides,x_ref[irefatom])

  # Move molecule to new position
  move!(x_new,newcm,beta,gamma,theta)

  nothing
end

