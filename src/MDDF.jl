module MDDF

  using Printf
  using Parameters
  using ProgressMeter
  using Statistics
  using FortranFiles
  using PDBTools

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
  include("./Selection.jl")
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
  # Default reading with the Chemfiles infrastructure
  Trajectory( filename :: String, solute :: Selection, solvent :: Selection, format="") =
    Chemfile(filename,solute,solvent,format)


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

  # Structures to report results
  include("./FileOperations.jl")
  include("./Result.jl")
  include("./Samples.jl")
  include("./finalresults.jl")
  include("./merge.jl")
  include("./save.jl")
  include("./read.jl")
  include("./write.jl")

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
  include("./mddf_linkedcells_self.jl")

  # for the linked-cell-parallel version
  include("./FrameData.jl")
  include("./sum.jl")
  include("./mddf_frame.jl")
  include("./mddf_frame_self.jl")
  include("./mddf_linkedcells_parallel.jl")

  # Parser to the default mddf method
  mddf = mddf_linkedcells_parallel

end




