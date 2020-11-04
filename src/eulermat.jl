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

  c1 = cos(beta) 
  s1 = sin(beta) 
  c2 = cos(gamma) 
  s2 = sin(gamma) 
  c3 = cos(theta) 
  s3 = sin(theta)

  A = SMatrix{3,3,Float64,9}( c2*c3,
                              c1*s3 + c3*s1*s2,
                              s1*s3 - c1*c3*s2,
                              -c2*s3,
                              c1*c3 - s1*s2*s3,
                              c1*s2*s3 + c3*s1,
                              s2,
                              -c2*s1,
                              c1*c2 )
  return A
end

