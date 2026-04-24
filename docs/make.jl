using Documenter
using ComplexMixtures
using PDBTools
using Plots
makedocs(
    modules=[
        ComplexMixtures,
        isdefined(Base, :get_extension) ? Base.get_extension(ComplexMixtures, :Plotting) : ComplexMixtures.Plotting,
    ],
    format = Documenter.HTML(top_menu = true),
    sitename="ComplexMixtures.jl",
    pages = [
        "Getting started" => Any[ 
            "Introduction" => "index.md",
            "Concepts" => "concepts.md",
            "Installation" => "installation.md",
            "Parallel execution" => "parallel.md",
            "Quick example" =>  "quickguide.md",
            "From Python" => "python.md",
        ],
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
        "References" => Any[
            "Primary citations" => "citations.md",
            "Applications" => "applications.md",
            "See also" => "seealso.md",
            "Breaking changes" => "breaking_changes.md",
        ],
    ],
)
deploydocs(
    repo="github.com/m3g/ComplexMixtures.jl.git",
    target="build",
    branch="gh-pages",
    versions=["stable" => "v^", "v#.#"],
)
