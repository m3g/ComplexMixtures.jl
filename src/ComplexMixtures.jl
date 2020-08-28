module ComplexMixtures

  using Printf
  using Parameters
  using ProgressMeter
  using Statistics
  using FortranFiles
  using PDBTools
  using StructTypes
  using JSON3
  using ThreadPools
  import Random

  # Input and Output data structures
  include("./FileOperations.jl")
  include("./OutputFiles.jl")
  include("./Density.jl")
  include("./Volume.jl")
  include("./Options.jl")

  # Structures used to store results
  include("./Result.jl")
  include("./Samples.jl")

  # Function to rigid-body move molecules
  include("./random.jl")
  include("./MoveAux.jl")
  include("./centerofcoordinates.jl")
  include("./eulermat.jl")
  include("./random_move.jl")
  include("./move.jl")
  include("./center_to_origin.jl")
  include("./writexyz.jl")

  # Structures and functions to deal with the solute and solvent selections
  include("./Selection.jl")
  include("./itype.jl")

  # Select solute and solvent using VMD
  include("./VMDselect.jl")

  # Function to print some data about the run
  include("./title.jl")

  # Structures and functions to read different types of trajectories
  include("./trajectory_formats/ChemFiles.jl")
  include("./trajectory_formats/NamdDCD.jl")
  include("./trajectory_formats/PDBTraj.jl")
  # Default reading with the Chemfiles infrastructure
  include("./Trajectory.jl")

  # Structure used for periodic boundary conditions
  include("./wrap.jl")
  include("./Box.jl")

  # Functions to compute distances
  include("./distance.jl")
  include("./minimumdistance.jl")

  # Functions to construct histograms
  include("./sphericalshellvolume.jl")
  include("./shellradius.jl")
  include("./sphereradiusfromshellvolume.jl")
  include("./setbin.jl")
  include("./viewmol.jl")
  include("./inbulk.jl")

  # Structures to report results
  include("./finalresults.jl")
  include("./merge.jl")
  include("./save.jl")
  include("./load.jl")
  include("./write.jl")
  include("./which_types.jl")
  include("./contrib.jl")

  # Implementation of mddf using naive algorithms
  include("./update_counters_frame.jl")
  include("./mddf_naive.jl")
  include("./mddf_naive_self.jl")
 
  # Structures and functions for the linked cell method
  include("./LinkedCells.jl")
  include("./DminMol.jl")
  include("./CutoffDistances.jl")
  include("./partialsort_cutoff.jl")
  include("./updatecounters.jl")
  include("./icell1D.jl")
  include("./wrap_cell.jl")
  include("./icell3D.jl")
  include("./initcells.jl")
  include("./increase_size.jl")
  include("./cutoffdcell.jl")
  include("./cutoffdistances.jl")
  include("./mddf_linkedcells.jl")
  include("./cutoffdcell_self.jl")
  include("./cutoffdistances_self.jl")

  # for the linked-cell-parallel version
  include("./FrameData.jl")
  include("./sum.jl")
  include("./mddf_frame.jl")
  include("./mddf_frame_self.jl")
  include("./mddf_linkedcells_parallel.jl")

  # Parser to the default mddf method
  include("./mddf_choose.jl")

end




