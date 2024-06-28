module Plotting

import Plots
import ComplexMixtures
using ComplexMixtures: Result, SoluteGroup, SolventGroup, contributions, ResidueContributions
using TestItems: @testitem
using PDBTools: Residue, residue_ticks, Atom, eachresidue, resnum

"""
    contourf(
        rc::ResidueContributions; 
        oneletter=false,
        xlabel="Residue",
        ylabel="r / Å",
        clims=nothing,
    )

Plot the contribution of each residue to the solute-solvent pair distribution function as a contour plot.
This function requires loading the `Plots` package.

# Arguments

- `rc::ResidueContributions`: The residue contributions to the solute-solvent pair distribution function,
   as computed by the `ResidueContributions` function.

# Optional arguments

- `oneletter::Bool`: Use one-letter residue codes. Default is `false`. One-letter codes are only available for the 20 standard amino acids.
- `xlabel` and `ylabel`: Labels for the x and y axes. Default is `"Residue"` and `"r / Å"`.
- `clims`: The color limits for the contour plot.

# Example

```julia-repl
julia> using ComplexMixtures, Plots, PDBTools

julia> results = load("mddf.json")

julia> atoms = readPDB("system.pdb", "protein")

julia> rc = ResidueContributions(results, atoms; oneletter=true)

julia> plt = contourf(rc)
```

This will produce a plot with the contribution of each residue to the solute-solvent pair distribution function,
as a contour plot, with the residues in the x-axis and the distance in the y-axis.

The resulting plot can customized using the standard mutating `plot!` function, for example, 

```julia-repl
julia> plot!(plt, size=(800, 400), title="Contribution per residue")
```

!!! compat
    This function requires loading the `Plots` package and is available in
    ComplexMixtures v2.5.0 or greater.

"""
function Plots.contourf(
    rc::ResidueContributions;
    oneletter=false,
    xlabel="Residue",
    ylabel="r / Å",
    clims=nothing,
    colorscale=:tempo,
)

    # Plot a contour courves with the density at each distance from each residue
    # colors, linewidths, etc. are defined here and can be tuned
    Plots.default(fontfamily="Computer Modern")
    if oneletter 
        xticks = (xticks[1], PDBTools.oneletter.(xticks[2][1:3]*xticks[2][4:end]))
    end
    plt = Plots.contourf(rc.xticks[1], rc.d, rc.residue_contributions,
        color=Plots.cgrad(colorscale),
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

function ComplexMixtures.contourf_per_residue(
    rc::ResidueContributions;
    oneletter=false,
    xlabel="Residue",
    ylabel="r / Å",
    clims=nothing,
)

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
    tmpplot = tempname() * ".png"
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
    tmpplot = tempname() * ".png"
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