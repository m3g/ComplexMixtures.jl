module ComplexMixtures

using PrecompileTools
using Printf
using ProgressMeter
using Statistics
using LinearAlgebra: norm, cross, dot, diag
using FortranFiles
using PDBTools
using StructTypes
using JSON3
using StaticArrays
using DocStringExtensions
using TestItems

import ChunkSplitters
import CellListMap
using .CellListMap.PeriodicSystems
import .CellListMap.PeriodicSystems: AbstractPeriodicSystem

import Random

export Selection, Trajectory, mddf, save, load, write, Options, Result
export merge
export overview
export VMDselect

function __init__()
    if !haskey(ENV, "COMPLEXMIXTURES_WARNING_v1.4") || ENV["COMPLEXMIXTURES_WARNING_v1.4"] == true
        @warn """\n
    
            ComplexMixtures - v1.4.2
    
            This version of ComplexMixtures is **not** the latest one. 
    
            If you did not intend to use this specific version, please update to the latest one by running:
    
            julia> import Pkg; Pkg.update("ComplexMixtures")
    
            One reason for incidentally installing this version is using a outdated version of Julia. 
            To avoid that, install Julia from the official website: https://julialang.org/downloads/
    
            If you intended to use this version, you can ignore this message, but follow the user guide
            with the correct version number, at: https://m3g.github.io/ComplexMixtures.jl/v1.4/ 

            To suppress this warning, use:

            julia> ENV["COMPLEXMIXTURES_WARNING_v1.4"] = false; using ComplexMixtures
            
        """ _file=nothing _line=nothing
    end
end


# Tools
export contributions, coordination_number, gr, grid3D

# Message for internal doc strings
const INTERNAL = "**Internal structure or function, interface may change.**"

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

# Functions to construct histograms
include("./viewmol.jl")

# Functions to compute distances
include("./minimum_distances.jl")

# Update counters routines
include("./updatecounters.jl")

# Main function
include("./mddf.jl")

# Comparison operators for ComplexMixtures types
include("./compare.jl")

# Tools
include("./tools/contributions.jl")
include("./tools/coordination_number.jl")
include("./tools/gr.jl")
include("./tools/grid3D.jl")

# Legacy interfaces
include("./legacy/legacy.jl")

# Precompilation directives
include("precompile.jl")

end # module ComplexMixtures
