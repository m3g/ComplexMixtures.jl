module MDDF

  using Printf
  using OffsetArrays

  include("./structures.jl")

  # Structures and functions to read different types of trajectories
  include("./trajectory_formats/NamdDCD.jl")
  include("./trajectory_formats/PDBTraj.jl")

  include("./Solute.jl")
  include("./Solvent.jl")

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

  include("./mddf.jl")

end




