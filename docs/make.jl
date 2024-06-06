import Pkg
Pkg.add("Documenter")
using Documenter
using ComplexMixtures
using PDBTools
using Plots
push!(LOAD_PATH, "../src/")
makedocs(
    modules = [
        ComplexMixtures, 
        isdefined(Base, :get_extension) ? Base.get_extension(ComplexMixtures, :Plotting) : ComplexMixtures.Plotting,
    ],
    sitename = "ComplexMixtures.jl",
    pages = [
        "Introduction" => "index.md",
        "Installation" => "installation.md",
        "Parallel execution" => "parallel.md",
        "Quick Guide" => "quickguide.md",
        "Examples:" => "examples.md",
        " ◦ Protein in water/glycerol" => "example1.md",
        " ◦ Polyacrylamide in DMF" => "example2.md",
        " ◦ POPC membrane in water/ethanol" => "example3.md",
        " ◦ Water/Glycerol mixture" => "example4.md",
        "Set solute and solvent" => "selection.md",
        "Loading the trajectory" => "trajectory.md",
        "Computing the MDDF" => "mddf.md",
        "Results" => "results.md",
        "Atomic and group contributions" => "contributions.md",
        "Save and load" => "save.md",
        "Multiple trajectories" => "multiple.md",
        "Options" => "options.md",
        "Tools" => "tools.md",
        "From Python" => "python.md",
        "Updating scripts" => "updating_scripts.md",
        "References" => "references.md",
    ],
)
deploydocs(
    repo = "github.com/m3g/ComplexMixtures.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
)
