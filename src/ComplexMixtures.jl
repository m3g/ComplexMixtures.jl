module ComplexMixtures

using Printf
using Parameters
using ProgressMeter
using Statistics
using FortranFiles
using PDBTools
using StructTypes
using JSON3
using StaticArrays
using DocStringExtensions
using TestItems

import CellListMap
using CellListMap.PeriodicSystems

import Random

export Selection, Trajectory, mddf, save, load, write, Options, Result
export contrib, merge
export overview, gr, grid3D
export VMDselect

# Message for internal doc strings
const INTERNAL = "Internal sructure or function, interface may change."

# Module for testing
include("./Testing.jl")

# Input and Output data structures
include("./io.jl")


# Options
include("./Options.jl")

# Structures and functions to deal with the solute and solvent selections
include("./VMDselect.jl")
include("./Selection.jl")

# Structures and functions to read different types of trajectories
include("./Trajectory.jl")

# Some functions to deal with rigid body calculations
include("./rigid_body.jl")

# Structures and functions to store and report results
include("./results.jl")

#= 
# Function to rigid-body move molecules
include("./random.jl")

# Functions to compute distances
include("./minimum_distances.jl")

# Functions to construct histograms
include("./viewmol.jl")

# Tools
include("./tools/gr.jl")
include("./tools/overview.jl")
include("./tools/grid3D.jl")
include("./isapprox.jl")

# Structures and functions for the linked cell method
include("./partialsort_cutoff.jl")
include("./updatecounters.jl")
include("./update_counters_frame.jl")

# for the linked-cell-parallel version
include("./FrameData.jl")
include("./mddf_frame.jl")
include("./mddf_frame_self.jl")
include("./mddf_linkedcells_parallel.jl")

# Main function
include("./mddf.jl")

=#

end
