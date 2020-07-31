module MDDF

  using Printf
  using Parameters
  using ProgressMeter
  using Statistics
  using FortranFiles

  include("./version.jl")

  # Function to rigid-body move molecules
  include("./MoveAux.jl")
  include("./centerofcoordinates.jl")
  include("./eulermat.jl")
  include("./random_move.jl")
  include("./move.jl")
  include("./center_to_origin.jl")
  include("./writexyz.jl")

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
  include("./trajectory_formats/ChemFiles.jl")
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
  include("./viewmol.jl")

  # Structures to report results
  include("./FileOperations.jl")
  include("./Result.jl")
  include("./Samples.jl")
  include("./finalresults.jl")

  # Implementation of mddf using naive algorithms
  include("./update_counters_frame.jl")
  include("./mddf_naive.jl")
  include("./mddf_naive_self.jl")
 
  # Structures and functions for the linked cell method
  dir="./linkedcells"
  include("$dir/LinkedCells.jl")
  include("$dir/DminMol.jl")
  include("$dir/CutoffDistances.jl")
  include("$dir/partialsort_cutoff.jl")
  include("$dir/updatecounters.jl")
  include("$dir/icell1D.jl")
  include("$dir/wrap_cell.jl")
  include("$dir/icell3D.jl")
  include("$dir/initcells.jl")
  include("$dir/increase_size.jl")
  include("$dir/cutoffdcell.jl")
  include("$dir/cutoffdistances.jl")
  include("$dir/mddf_linkedcells.jl")
  include("$dir/cutoffdcell_self.jl")
  include("$dir/cutoffdistances_self.jl")
  include("$dir/mddf_linkedcells_self.jl")

  # for the linked-cell-parallel version
  dir="./linkedcells_parallel"
  include("$dir/FrameData.jl")
  include("$dir/sum.jl")
  include("$dir/mddf_frame.jl")
  include("$dir/mddf_frame_self.jl")
  include("$dir/mddf_linkedcells_parallel.jl")

  # Parser to the default mddf method
  mddf = mddf_linkedcells_parallel

end




