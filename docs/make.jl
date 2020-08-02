using Documenter, MDDF
push!(LOAD_PATH,"../src/")
makedocs(
    modules=[MDDF],
    sitename="MDDF.jl",
    pages = [
        "Introduction" => "index.md",
        "Installation" => "installation.md",
        "Quick Guide" => "quickguide.md",
        "Set solute and solvent" => "selection.md",
        "Loading the trajectory" => "trajectory.md",
        "Computing the MDDF" => "mddf.md",
        "Results" => "results.md",
        "Multiple trajectories" => "multiple.md",
        "Parallel execution" => "parallel.md",
        "Options" => "options.md",
        "Examples" => "examples.md",
        "References" => "references.md"
    ]
)



