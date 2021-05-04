"""

Structure to contain data needed to compute the mddf for a single frame

"""
mutable struct FrameData{T<:Trajectory,V}

  trajectory::T
  volume_frame::Volume
  rdf_count_random_frame::Vector{Float64}
  md_count_random_frame::Vector{Float64}
  
  dc::CutoffDistances
  dmin_mol::Vector{DminMol}
  dref_mol::Vector{Float64}
  x_solvent_random::Vector{V}

  lc_solvent::LinkedCells

end

function FrameData(trajectory::Trajectory, R::Result ) 
  return FrameData(trajectory,                                       # trajectory
                   Volume(R.nbins),                                  # volume_frame
                   zeros(R.nbins),                                   # rdf_count_random_frame
                   zeros(R.nbins),                                   # md_count_random_frame
                   CutoffDistances(trajectory.solvent.natoms),       # dc       
                   [ DminMol(+Inf,i,0,0) for i in 1:trajectory.solvent.nmols ], # dmin_mol
                   zeros(trajectory.solvent.nmols),                             # dref_mol
                   similar(trajectory.x_solvent),                               # x_solvent_random
                   LinkedCells(trajectory.solvent.natoms))                      # lc_solvent
end

