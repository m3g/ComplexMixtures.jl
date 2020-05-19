#
# Structure used for linked cells in computation of small distances
#

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


