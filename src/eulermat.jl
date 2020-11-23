#
# function eulermat: This routine was added because it defines 
#                    the rotation in the "human" way, an is thus used
#                    to set the position of the fixed molecules. 
#     That means: beta is a counterclockwise rotation around x axis.
#                 gamma is a counterclockwise rotation around y axis.
#                 theta is a counterclockwise rotation around z axis.
# 

function eulermat( beta, gamma, theta, deg :: String )

   if deg != "degree" 
     error("ERROR: to use radians just omit the last parameter")
   end

   beta = beta*π/180
   gamma = gamma*π/180
   theta = theta*π/180
   A = eulermat( beta, gamma, theta )
   return A

end

function eulermat(beta :: T, gamma :: T, theta ::T) where T
  c1 = cos(beta) 
  s1 = sin(beta) 
  c2 = cos(gamma) 
  s2 = sin(gamma) 
  c3 = cos(theta) 
  s3 = sin(theta)
  A = SMatrix{3,3,T,9}( c2*c3,
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

