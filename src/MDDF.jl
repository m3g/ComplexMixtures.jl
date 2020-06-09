module MDDF

  using Printf
  using OffsetArrays

  include("./setbin.jl")

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
  include("./LinkedLists.jl")

  # Structures and functions to read different types of trajectories
  include("./trajectory_formats/NamdDCD.jl")
  include("./trajectory_formats/PDBTraj.jl")

  include("./VMDselect.jl")

  include("./FileOperations.jl")
  include("./format.jl")

  include("./compcart.jl")
  include("./centerofcoordinates.jl")
  include("./eulermat.jl")

  include("./phantomcells.jl")
  include("./movephantomcoor.jl")
  include("./pbc.jl")
  include("./smalldcell.jl")
  include("./smalldistances.jl")

  include("./sphericalshellvolume.jl")
  include("./shellradius.jl")
  include("./sphereradiusfromshellvolume.jl")


  #include("./mddf.jl")

  include("./mddf_naive.jl")
  include("./mddf_naive_single.jl")
  include("./mddf_naive_homogeneous.jl")

end




