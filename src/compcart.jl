#
# function compcart: Compute cartesian coordinates using
#                    the center of mass, the canonical coordinates
#                    and the rotation matrix
#
# xcm: Coordinates of the center of mass of the molecule
# A : rotation matrix
# xref: reference coordinates (assumed to be centered at the origin)
#      

# From the rotation matrix

function compcart( xcm :: Vector{Float64}, 
                   A :: Matrix{Float64}, 
                   xref :: Vector{Float64} )
  x = Vector{Float64}(undef,3)
  x[1] = xcm[1] + xref[1]*A[1,1] + xref[2]*A[2,1] + xref[3]*A[3,1]    
  x[2] = xcm[2] + xref[1]*A[1,2] + xref[2]*A[2,2] + xref[3]*A[3,2]    
  x[3] = xcm[3] + xref[1]*A[1,3] + xref[2]*A[2,3] + xref[3]*A[3,3]    
  return x
end

# From the angles

function compcart( xcm :: Vector{Float64}, 
                   xref :: Vector{Float64}, 
                   beta :: Float64, 
                   gamma :: Float64, 
                   theta :: Float64 ) 
   A = eulermat(beta, gamma, theta)
   x = compcart( xcm, A, xref )
   return x
end

# From the angles, but angles given in degrees instead of radians

function compcart( xcm :: Vector{Float64}, 
                   xref :: Vector{Float64}, 
                   beta :: Float64, 
                   gamma :: Float64, 
                   theta :: Float64, 
                   deg :: String ) 
   A = eulermat(beta, gamma, theta, deg)
   x = compcart( xcm, A, xref )
   return x
end

# For a set of atoms instead a single atom, modifies xnew

function compcart!( natoms,
                    xcm :: Vector{Float64}, 
                    xref :: Vector{Float64}, 
                    yref :: Vector{Float64}, 
                    zref :: Vector{Float64}, 
                    beta :: Float64, 
                    gamma :: Float64, 
                    theta :: Float64, 
                    xnew :: Vector{Float64},  
                    ynew :: Vector{Float64},  
                    znew :: Vector{Float64} ) 
   A = eulermat(beta, gamma, theta)
   for i in 1:natoms
     xin = [ xref[i], yref[i], zref[i] ]
     x = compcart( xcm, A, xin )
     xnew[i] = x[1]
     ynew[i] = x[2]
     znew[i] = x[3]
   end
end



