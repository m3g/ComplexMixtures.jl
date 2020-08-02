using Documenter, MDDF
push!(LOAD_PATH,"../src/")
makedocs(
    modules= [MDDF],
    sitename="MDDF.jl",
    pages = [
        "Introduction" => "index.md",
        "Installation" => "installation.md",
        "Quick Guide" => "quickguide.md",
        "Selections" => "selection.md",
        "Multiple trajectories" => "multiple.md",
        "Save results" => "save.md",
        "Parallel execution" => "parallel.md",
        "Options" => "options.md",
        "References" => "references.md"
    ]
)



