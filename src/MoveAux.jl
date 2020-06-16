#
# Auxiliary structure for random rotations and translations
#

struct MoveAux
  oldcm :: Vector{Float64}
  newcm :: Vector{Float64}
  angles :: Vector{Float64}
  A :: Matrix{Float64}
  x :: Matrix{Float64}
end
MoveAux(n) = MoveAux(zeros(3), zeros(3), zeros(3), zeros(3,3), zeros(n,3))

