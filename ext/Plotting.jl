module Plotting

using TestItems: @testitem
import Plots
import ComplexMixtures
import PDBTools
using ComplexMixtures: Result, SoluteGroup, SolventGroup, contributions, ResidueContributions
using PDBTools: Residue, residue_ticks, Atom, eachresidue, resnum

"""
    contourf(
        rc::ResidueContributions; 
        oneletter=false,
        xlabel="Residue",
        ylabel="r / Å",
        clims=nothing,
        colorscale=:tempo,
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
- `colorscale`: The color scale for the contour plot. Default is `:tempo`. We suggest `:bwr` if the zero
  is at the middle of the color scale (use `clims` to adjust the color limits).  

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
    xticks_range=nothing,
    oneletter=false,
    xlabel="Residue",
    ylabel="r / Å",
    clims=nothing,
    colorscale=:tempo,
)

    residue_numbers = [ parse(Int, label[4:end]) for label in rc.xticks[2] ]

    # Check if the range is fine
    xticks_range, tick_range = if isnothing(xticks_range)
        first(rc.xticks[1]):1:last(rc.xticks[1]), 1:1:length(rc.xticks[1])
    else
        ifirst = findfirst(==(xticks_range[1]), rc.xticks[1])
        ilast = findfirst(==(xticks_range[end]), rc.xticks[1])
        if isnothing(ifirst) || isnothing(ilast)
            throw(ArgumentError("""\n
    
                The xticks_range $xticks_range is out of bounds.
    
                The first and last residue numbers of the selection are: $(first(residue_numbers)) and $(last(residue_numbers)).
                
            """
            ))
        end
        xticks_range, ifirst:step(xticks_range):ilast
    end

    # Plot a contour courves with the density at each distance from each residue
    # colors, linewidths, etc. are defined here and can be tuned
    tick_marks = rc.xticks[1][tick_range]
    tick_labels = rc.xticks[2][tick_range]
    xticks = if oneletter 
        for i in eachindex(tick_labels)
            tick_labels[i] = PDBTools.oneletter(tick_labels[i][1:3]) * tick_labels[i][4:end] 
        end 
        (tick_marks, tick_labels) 
    else
        (tick_marks, tick_labels) 
    end
    Plots.default(fontfamily="Computer Modern")
    plt = Plots.contourf(residue_numbers, rc.d, rc.residue_contributions,
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

ComplexMixtures.contourf_per_residue(rc::ResidueContributions; kwargs...) = Plots.contourf(rc; kwargs...)

"""
    contourf_per_residue(
        results::Result, atoms::AbstractVector{PDBTools.Atom}; 
        xticks_range=nothing, 
        dmin=1.5, dmax=3.5, 
        oneletter=false,
        xlabel="Residue",
        ylabel="r / Å",
        type=:mddf,
        clims=nothing,
        colorscale=:tempo,
    )

Plot the contribution of each residue to the solute-solvent pair distribution function as a contour plot.
This function requires loading the `Plots` package.

# Arguments

- `results::Result`: The result of a `mddf` call.
- `atoms::AbstractVector{Atom}`: The atoms of the solute.

# Optional arguments

- `xticks_range`: The range of residues to plot. Default is all residues. Use
  a step to plot every `n` residues, e.g., `1:2:30`.
- `dmin::Real`: The minimum distance to plot. Default is `1.5`.
- `dmax::Real`: The maximum distance to plot. Default is `3.5`.
- `oneletter::Bool`: Use one-letter residue codes. Default is `false`. One-letter codes are only available for the 20 standard amino acids.
- `xlabel` and `ylabel`: Labels for the x and y axes. Default is `"Residue"` and `"r / Å"`.
- `type::Symbol`: That data to plot. Default is `:mddf` for MDDF contributions. Options are `:coordination_number`, and `:mddf_count`.
- `clims`: The color limits for the contour plot.
- `colorscale`: The color scale for the contour plot. Default is `:tempo`. We suggest `:bwr` if the zero
  is at the middle of the color scale (use `clims` to adjust the color limits).

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

    The `type` and `clims` arguments, and the support for a step in `xticks_range` were introduced in ComplexMixtures v2.3.0.

"""
function ComplexMixtures.contourf_per_residue(
    results::Result, atoms::AbstractVector{Atom}; 
    xticks_range=nothing,
    dmin=1.5, dmax=3.5,
    type=:mddf,
    kargs...
)
    first_atom = findfirst(at -> resnum(at) == first(xticks_range), atoms)
    last_atom = findfirst(at -> resnum(at) == last(xticks_range), atoms)
    if isnothing(first_atom) || isnothing(last_atom)
        throw(ArgumentError("""\n

            The xticks_range $xticks_range is out of bounds.

            The first and last residue numbers of the selection are: $(resnum(first(atoms))) and $(resnum(last(atoms))).

        """
        ))
    end
    rc = ResidueContributions(results, atoms[first_atom:last_atom]; dmin, dmax, type)
    return Plots.contourf(rc; xticks_range, kargs...)
end

@testitem "contourf/contourf_per_residue" begin
    using ComplexMixtures
    using ComplexMixtures.Testing
    using Plots
    using PDBTools
    # Load example output file (computed in the previous script)
    protein = readPDB(joinpath(Testing.data_dir, "Gromacs/system.pdb"), "protein")
    results = load(joinpath(Testing.data_dir, "Gromacs/protein_EMI.json"))

    rc = ResidueContributions(results, select(protein, "residue >= 50 and residue <= 75"))
    plt = contourf(rc)
    tmpplot = tempname() * ".png"
    savefig(plt, tmpplot)
    @test isfile(tmpplot)

    plt = contourf_per_residue(results, protein; xticks_range=50:75, oneletter=true)
    tmpplot = tempname() * ".png"
    savefig(plt, tmpplot)
    @test isfile(tmpplot)

    plt = contourf_per_residue(results, protein; xticks_range=50:75)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)

    plt = contourf_per_residue(results, protein; xticks_range=50:2:70)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)

    @test_throws ArgumentError contourf_per_residue(results, SoluteGroup("ALA"))
    @test_throws ArgumentError contourf_per_residue(results, protein; xticks_range=50:80)
end

end # module Plotting