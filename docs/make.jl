#using DocumenterLaTeX
using ComplexMixtures
using Documenter
using ComplexMixtures
push!(LOAD_PATH,"../src/")
makedocs(
    modules=[ComplexMixtures],
    sitename="ComplexMixtures",
    pages = [
        "Introduction" => "index.md",
        "Installation" => "installation.md",
        "Quick Guide" => "quickguide.md",
        "Set solute and solvent" => "selection.md",
        "Loading the trajectory" => "trajectory.md",
        "Computing the MDDF" => "mddf.md",
        "Results" => "results.md",
        "Atomic and group contributions" => "contrib.md",
        "Save and load" => "save.md",
        "Multiple trajectories" => "multiple.md",
        "Parallel execution" => "parallel.md",
        "Options" => "options.md",
        "Tools" => "tools.md",
        "Examples" => "examples.md",
        "References" => "references.md"
    ]
)



