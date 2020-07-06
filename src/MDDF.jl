module MDDF

  using Printf
  using OffsetArrays
  using Parameters
  using Statistics
  using ProgressMeter

  include("./version.jl")

  # Function to rigid-body move molecules
  include("./MoveAux.jl")
  include("./centerofcoordinates.jl")
  include("./eulermat.jl")
  include("./random_move.jl")
  include("./move.jl")

  # Structures and functions to deal with the solute and solvent selections
  include("./SoluteOrSolvent.jl")
  include("./itype.jl")

  # Select solute and solvent using VMD
  include("./VMDselect.jl")

  # Input and Output data structures
  include("./OutputFiles.jl")
  include("./Density.jl")
  include("./Volume.jl")
  include("./Options.jl")

  # Structures and functions to read different types of trajectories
  include("./trajectory_formats/NamdDCD.jl")
  include("./trajectory_formats/PDBTraj.jl")

  # Wrapping functions
  include("./wrap.jl")

  # Structure used for periodic boundary conditions
  include("./linkedcells/Box.jl")

  # Functions to compute distances
  include("./distance.jl")
  include("./minimumdistance.jl")

  # Functions to construct histograms
  include("./sphericalshellvolume.jl")
  include("./shellradius.jl")
  include("./sphereradiusfromshellvolume.jl")
  include("./setbin.jl")

  # Structures to report results
  include("./FileOperations.jl")
  include("./Result.jl")
  include("./finalresults.jl")
  include("./finalresults_self.jl")

  # Implementation of mddf using naive algorithms
  include("./mddf_naive.jl")
  include("./mddf_naive_self.jl")
 
  # Structures and functions for the linked cell method

  #include("./linkedcells_v1/LinkedCells.jl")
  #include("./linkedcells_v1/icell3D.jl")
  #include("./linkedcells_v1/ijkcell.jl")
  #include("./linkedcells_v1/initcells.jl")
  #include("./linkedcells_v1/CutoffDistances.jl")
  #include("./linkedcells_v1/cutoffdcell.jl")
  #include("./linkedcells_v1/cutoffdistances.jl")
  #include("./linkedcells_v1/keepunique.jl")
  #include("./linkedcells_v1/keepminimum.jl")
  #include("./linkedcells_v1/mddf_linkedcells.jl")

  include("./linkedcells/LinkedCells.jl")
  include("./linkedcells/DminMol.jl")
  include("./linkedcells/partialsort_cutoff.jl")
  include("./linkedcells/updatecounters.jl")
  include("./linkedcells/icell1D.jl")
  include("./linkedcells/wrap_cell.jl")
  include("./linkedcells/icell3D.jl")
  include("./linkedcells/initcells.jl")
  include("./linkedcells/CutoffDistances.jl")
  include("./linkedcells/cutoffdcell.jl")
  include("./linkedcells/cutoffdistances.jl")
  include("./linkedcells/keepunique.jl")
  include("./linkedcells/keepminimum.jl")
  include("./linkedcells/mddf_linkedcells.jl")

  # Parser to the default mddf method in each case
  #include("./mddf.jl")

end




