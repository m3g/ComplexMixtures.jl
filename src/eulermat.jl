#
# function eulermat: This routine was added because it defines 
#                    the rotation in the "human" way, an is thus used
#                    to set the position of the fixed molecules. 
#     That means: beta is a counterclockwise rotation around x axis.
#                 gamma is a counterclockwise rotation around y axis.
#                 theta is a counterclockwise rotation around z axis.
# 

function eulermat( beta :: Float64,
                   gamma :: Float64,
                   theta :: Float64, 
                   deg :: String )

   if deg != "degree" 
     error("ERROR: to use radians just omit the last parameter")
   end

   beta = pi * beta / 180.
   gamma = pi * gamma / 180.
   theta = pi * theta / 180.
   A = eulermat( beta, gamma, theta )
   return A

end

function eulermat( beta :: Float64,
                   gamma :: Float64,
                   theta :: Float64 )

  A = Matrix{Float64}(undef,3,3)

  c1 = cos(beta) 
  s1 = sin(beta) 
  c2 = cos(gamma) 
  s2 = sin(gamma) 
  c3 = cos(theta) 
  s3 = sin(theta)

  A[1,1] = c2*c3
  A[1,2] = c1*s3 + c3*s1*s2
  A[1,3] = s1*s3 - c1*c3*s2

  A[2,1] = -c2*s3
  A[2,2] = c1*c3 - s1*s2*s3
  A[2,3] = c1*s2*s3 + c3*s1

  A[3,1] = s2
  A[3,2] = -c2*s1
  A[3,3] = c1*c2         

  return A

end

# Function that performs the same computation, but uptating the provided
# auxiliary structure, to avoid new allocations

function eulermat!( aux :: MoveAux )

  c1 = cos(aux.angles[1]) 
  s1 = sin(aux.angles[1]) 
  c2 = cos(aux.angles[2]) 
  s2 = sin(aux.angles[2]) 
  c3 = cos(aux.angles[3]) 
  s3 = sin(aux.angles[3])

  aux.A[1,1] = c2*c3
  aux.A[1,2] = c1*s3 + c3*s1*s2
  aux.A[1,3] = s1*s3 - c1*c3*s2

  aux.A[2,1] = -c2*s3
  aux.A[2,2] = c1*c3 - s1*s2*s3
  aux.A[2,3] = c1*s2*s3 + c3*s1

  aux.A[3,1] = s2
  aux.A[3,2] = -c2*s1
  aux.A[3,3] = c1*c2         

  return nothing
end
 
