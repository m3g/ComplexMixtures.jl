using Documenter, MDDF
push!(LOAD_PATH,"../src/")
makedocs(
    modules= [MDDF],
    sitename="MDDF.jl",
    pages = [
        "Introduction" => "index.md",
        "Selections" => "selection.md"
    ]

)



