#
# Structure that contains which are the atoms belonging to the solute and 
# the solvent (groups 1 and 2, respectively)
#

struct Groups

  n1 :: Int64
  n2 :: Int64
  group1 :: Vector{Int}
  group2 :: Vector{Int}
  group1_box :: Vector{Int64}
  group2_box :: Vector{Int64}

end

