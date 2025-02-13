import Pkg
Pkg.add("Documenter")
using Documenter
using ComplexMixtures
using PDBTools
using Plots
push!(LOAD_PATH, "../src/")
makedocs(
    modules=[
        ComplexMixtures,
        isdefined(Base, :get_extension) ? Base.get_extension(ComplexMixtures, :Plotting) : ComplexMixtures.Plotting,
    ],
    sitename="ComplexMixtures.jl",
    pages=[
        "Introduction" => "index.md",
        "Install and run" => Any[ 
            "Installation" => "installation.md",
            "Parallel execution" => "parallel.md",
            "From Python" => "python.md",
        ],
        "Quick Guide" => "quickguide.md",
        "Examples" => Any[
            "Running the examples" => "examples.md",
            "Protein in water/glycerol" => "example1.md",
            "Polyacrylamide in DMF" => "example2.md",
            "POPC membrane in water/ethanol" => "example3.md",
            "Water/Glycerol mixture" => "example4.md",
        ],
        "Setup and run" => Any[
            "Set solute and solvent" => "selection.md",
            "Computing the MDDF" => "mddf.md",
            "Save and load" => "save.md",
            "Multiple trajectories" => "multiple.md",
            "Options" => "options.md",
        ],
        "Analysis" => Any[
            "Results" => "results.md",
            "Atomic and group contributions" => "contributions.md",
            "Density maps" => "density_maps.md",
            "Coordination numbers" => "coordination_numbers.md",
            "Tools" => "tools.md",
        ],
#        "Updating scripts" => "updating_scripts.md",
        "References" => "references.md",
    ],
)
deploydocs(
    repo="github.com/m3g/ComplexMixtures.jl.git",
    target="build",
    branch="gh-pages",
    versions=["stable" => "v^", "v#.#"],
)
