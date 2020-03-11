#
# Structure used for linked cells in computation of small distances
#

struct Frame
  
  sides :: Vector{Float64}
  x :: Vector{Float64}
  y :: Vector{Float64}
  z :: Vector{Float64}

end

struct Groups

  n1 :: Int64
  n2 :: Int64
  group1 :: Vector{Int}
  group2 :: Vector{Int}
  group1_box :: Vector{Int64}
  group2_box :: Vector{Int64}

end

mutable struct SmalldLinkedLists

  iatomfirst :: Array{Int64}
  iatomnext :: Vector{Int64}
  nboxes :: Vector{Int64}
  dbox :: Vector{Float64}
  nbdim :: Vector{Int64}
  cutoff :: Float64

end

mutable struct SmallDistances

  n :: Int64
  nmax :: Int64
  index :: Vector{Int64}
  d :: Vector{Int64}
  memerror :: Bool

end

struct DistanceData

  frame :: Frame
  groups :: Groups
  lists :: SmalldLinkedLists
  smalld :: SmallDistances

end


