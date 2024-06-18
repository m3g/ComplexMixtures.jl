module Plotting

import Plots
import ComplexMixtures
using ComplexMixtures: Result, SoluteGroup, SolventGroup, contributions
using TestItems: @testitem
using PDBTools: Residue, residue_ticks, Atom, eachresidue, resnum

"""
    contourf_per_residue(
        results::Result, atoms::AbstractVector{PDBTools.Atom}; 
        residue_range=1:1:length(eachresidue(atoms)), 
        dmin=1.5, dmax=3.5, 
        oneletter=false,
        xlabel="Residue",
        ylabel="r / Å",
        type=:mddf,
        clims=nothing,
    )

Plot the contribution of each residue to the solute-solvent pair distribution function as a contour plot.
This function requires loading the `Plots` package.

# Arguments

- `results::Result`: The result of a `mddf` call.
- `atoms::AbstractVector{Atom}`: The atoms of the solute.

# Optional arguments

- `residue_range`: The range of residues to plot. Default is `1:1:length(eachresidue(atoms))`. Use
  a step to plot every `n` residues, e.g., `1:2:30`.
- `dmin::Real`: The minimum distance to plot. Default is `1.5`.
- `dmax::Real`: The maximum distance to plot. Default is `3.5`.
- `oneletter::Bool`: Use one-letter residue codes. Default is `false`. One-letter codes are only available for the 20 standard amino acids.
- `xlabel` and `ylabel`: Labels for the x and y axes. Default is `"Residue"` and `"r / Å"`.
- `type::Symbol`: That data to plot. Default is `:mddf` for MDDF contributions. Options are `:coordination_number`, and `:mddf_count`.
- `clims`: The color limits for the contour plot.

# Example

```julia-repl
julia> using ComplexMixtures, Plots, PDBTools

julia> results = load("mddf.json")

julia> atoms = readPDB("system.pdb", "protein")

julia> plt = contourf_per_residue(results, atoms; oneletter=true)
```

This will produce a plot with the contribution of each residue to the solute-solvent pair distribution function,
as a contour plot, with the residues in the x-axis and the distance in the y-axis.

The resulting plot can customized using the standard mutating `plot!` function, for example, 

```julia-repl
julia> plot!(plt, size=(800, 400), title="Contribution per residue")
```

!!! compat
    This function requires loading the `Plots` package and is available in 
    ComplexMixtures v2.2.0 or greater.

    The `type` and `clims` arguments, and the support for a step in `residue_range` were introduced in ComplexMixtures v2.3.0.

"""
function ComplexMixtures.contourf_per_residue(
    results::Result, atoms::AbstractVector{Atom};
    residue_range=1:length(eachresidue(atoms)),
    dmin=1.5, dmax=3.5,
    oneletter=false,
    xlabel="Residue",
    ylabel="r / Å",
    type=:mddf,
    clims=nothing,
)

    # collect the list of residues (using PDBTools)
    residues = collect(eachresidue(atoms))

    # Check if the range is fine
    if last(residue_range) > length(residues) || first(residue_range) < 1
        throw(ArgumentError("""\n
        
            The residue_range $residue_range is out of bounds

            The atom selection provided has $(length(residues)) residues. 
            Select a range where the first and last indices are within 1 and $(length(residues)).
            
        """
        ))
    end

    # Create matrix that will cotain the contribution per 
    # residue as a function of the distance:
    # number of rows of the mddf histogram is (length(results.d)) and 
    # number of columns equal to the number of residues
    rescontrib = zeros(length(results.d), length(residues))

    # Each column is then filled up with the contributions of each residue
    for (ires, residue) in enumerate(residues)
        rescontrib[:, ires] .= contributions(results, SoluteGroup(residue); type)
    end

    # Plot only for distances within 1.5 and 3.5:
    # find the indexes of the distances that are within the range
    idmin = findfirst(d -> d > dmin, results.d)
    idmax = findfirst(d -> d > dmax, results.d)

    # Obtain pretty labels for the residues in the x-axis (using PDBTools)
    xticks = residue_ticks(atoms; 
        first=residue_range[1],
        last=residue_range[end],
        stride=step(residue_range),
        oneletter
    )

    @show xticks
    # Plot a contour courves with the density at each distance from each residue
    # colors, linewidths, etc. are defined here and can be tuned
    Plots.default(fontfamily="Computer Modern")
    plt = Plots.contourf(residue_range, results.d[idmin:idmax], rescontrib[idmin:idmax, residue_range],
        color=Plots.cgrad(:tempo), 
        linewidth=1, 
        linecolor=:black,
        colorbar=:none, 
        levels=5,
        xlabel=xlabel, 
        ylabel=ylabel,
        xticks=xticks, xrotation=60, xtickfont=Plots.font(8, "Computer Modern"),
        size=(700, 400),
        margin=0.5Plots.PlotMeasures.cm,
        framestyle=:box,
        clims=clims,
    )

    # return the plot
    return plt
end

# Custom error message for common mistaken call
function ComplexMixtures.contourf_per_residue(result, g::Union{SoluteGroup,SolventGroup}, args...; kwargs...)
    throw(ArgumentError("""\n

        contourf_per_residue cannot be run if the MDDFs were computed with custom groups.

        This is because it requires the contribution of each independent atom to accumulate
        the residue contributions. The second argument of the function must be the vector of 
        atoms of the solute, not a group selection.

        Please read the documention, by typing: `? contourf_per_residue`
    
    """)) 
end

@testitem "contourf_per_residue" begin 
    using ComplexMixtures
    using ComplexMixtures.Testing
    using Plots
    using PDBTools
    # Load example output file (computed in the previous script)
    protein = readPDB(joinpath(Testing.data_dir, "Gromacs/system.pdb"), "protein")
    results = load(joinpath(Testing.data_dir, "Gromacs/protein_EMI.json"))
    plt = contourf_per_residue(results, protein; residue_range=50:75, oneletter=true)
    tmpplot = tempname()*".png"
    savefig(plt, tmpplot)
    @test isfile(tmpplot)
    plt = contourf_per_residue(results, protein; residue_range=50:75)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)
    plt = contourf_per_residue(results, protein; residue_range=50:2:70)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)
    @test_throws ArgumentError contourf_per_residue(results, SoluteGroup("ALA"))
    @test_throws ArgumentError contourf_per_residue(results, protein; residue_range=50:80)
end

end