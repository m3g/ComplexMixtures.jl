module MDDF

  using Printf
  using OffsetArrays
  using Parameters
  using Statistics
  using ProgressMeter

  include("./setbin.jl")
  include("./MoveAux.jl")
  include("./random_move.jl")
  include("./move.jl")

  # Input and Output data structures
  include("./OutputFiles.jl")
  include("./SoluteOrSolvent.jl")
  include("./Density.jl")
  include("./Volume.jl")
  include("./Options.jl")
  include("./Result.jl")

  # Structures to contain data critical for calculation performance
  include("./Frame.jl")
  include("./Groups.jl")

  # Structures and functions to read different types of trajectories
  include("./trajectory_formats/NamdDCD.jl")
  include("./trajectory_formats/PDBTraj.jl")

  include("./VMDselect.jl")

  include("./FileOperations.jl")
  include("./format.jl")

  include("./compcart.jl")
  include("./centerofcoordinates.jl")
  include("./eulermat.jl")

  include("./wrap.jl")

  include("./sphericalshellvolume.jl")
  include("./shellradius.jl")
  include("./sphereradiusfromshellvolume.jl")

  include("./distance.jl")
  include("./minimumdistance.jl")

  include("./finalresults.jl")
  include("./finalresults_self.jl")

  #include("./mddf.jl")
  include("./mddf_naive.jl")
  include("./mddf_naive_self.jl")
 
  # Structures and functions for the linked cell method
  include("./linkedcells/LinkedCells.jl")
  include("./linkedcells/icell3D.jl")
  include("./linkedcells/droprepeated.jl")
  include("./linkedcells/initcells.jl")
  include("./linkedcells/first_atom_in_cell.jl")

  #include("./mddf_lindedcells.jl")

end




