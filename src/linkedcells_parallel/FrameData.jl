#
# Structure to contain data needed to compute the mddf for a single frame
#

mutable struct FrameData

  trajectory # trajectory format
  volume_frame :: Volume
  rdf_count_random_frame :: Vector{Float64}
  sides :: Vector{Float64}
  box :: Box
  solute_center :: Vector{Float64}
  
  dc :: CutoffDistances
  dmin_mol :: Vector{DminMol}
  dref_mol :: Vector{Float64}
  x_solvent_random :: Array{Float64} 

  lc_solvent :: LinkedCells

  moveaux :: MoveAux
  nsamples :: Int64

end
