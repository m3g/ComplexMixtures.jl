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
using StaticArrays
import Random

export Selection, Trajectory, mddf, save, load, write, Options
export contrib, merge 
export overview, gr, grid3D
export VMDSelect

# Input and Output data structures
include("./FileOperations.jl")
include("./OutputFiles.jl")
include("./Density.jl")
include("./Volume.jl")
include("./Options.jl")
include("./Units.jl")

# Function to rigid-body move molecules
include("./random.jl")
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

# Structures and functions to read different types of trajectories
include("./Trajectory.jl")
include("./trajectory_formats/ChemFiles.jl")
include("./trajectory_formats/NamdDCD.jl")
include("./trajectory_formats/PDBTraj.jl")

# Structures used to store results
include("./isautocorrelation.jl")
include("./SolSummary.jl")
include("./Result.jl")
include("./Samples.jl")

# Function to print some data about the run
include("./title.jl")

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

# Tools
include("./tools/gr.jl")
include("./tools/overview.jl")
include("./tools/grid3D.jl")
include("./isapprox.jl")

# Structures and functions for the linked cell method
include("./LinkedCells.jl")
include("./DminMol.jl")
include("./CutoffDistances_Struct.jl")
include("./partialsort_cutoff.jl")
include("./updatecounters.jl")
include("./update_counters_frame.jl")
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

# Main function
include("./mddf.jl")

end




