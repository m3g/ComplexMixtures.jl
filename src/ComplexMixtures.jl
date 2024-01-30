module ComplexMixtures

import CellListMap
import .CellListMap.PeriodicSystems
import ChunkSplitters
import JSON3
import PDBTools
import PrecompileTools
import Random
import StructTypes

using .CellListMap.PeriodicSystems: AbstractPeriodicSystem, PeriodicSystem, 
    update_unitcell!, map_pairwise!
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using LinearAlgebra: norm, cross, dot, diag
using FortranFiles: FortranFile, rewind, FString, read, close
using Printf: @sprintf
using ProgressMeter: Progress, next!
using StaticArrays: SVector, SMatrix, @SMatrix, MMatrix
using Statistics: mean, std
using TestItems: @testitem


# Data types
export Trajectory, Options, Result 
export AtomSelection, SoluteGroup, SolventGroup

# Functions
export mddf
export overview
export save, load, write, merge
export atom_group, atom_group_name, atom_group_names

# Tools
export contributions, coordination_number, gr, grid3D

# Module for testing
include("./Testing.jl")

# Input and Output data structures
include("./io.jl")

# Options
include("./Options.jl")

# Structures and functions to deal with the solute and solvent selections
include("./AtomSelection.jl")

# Structures and functions to read different types of trajectories
include("./Trajectory.jl")

# Some functions to deal with rigid body calculations
include("./rigid_body.jl")

# Functions to construct histograms
include("./viewmol.jl")

# Structures and functions to store and report results
include("./results.jl")

# Comparison operators for ComplexMixtures types
include("./compare.jl")

# Functions to compute distances
include("./minimum_distances.jl")

# Update counters routines
include("./updatecounters.jl")

# Main function
include("./mddf.jl")

# Tools
include("./tools/contributions.jl")
include("./tools/coordination_number.jl")
include("./tools/gr.jl")
include("./tools/grid3D.jl")

# Legacy interfaces
include("./legacy/legacy.jl")

# Precompilation directives
#include("precompile.jl")

end # module ComplexMixtures
