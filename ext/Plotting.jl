module Plotting

using TestItems: @testitem
import Plots
import ComplexMixtures
import PDBTools
using ComplexMixtures: Result, SoluteGroup, SolventGroup, contributions,
    ResidueContributions, _set_clims_and_colorscale!
using PDBTools: Residue, residue_ticks, Atom, eachresidue, resnum

"""
    contourf(
        rc::ResidueContributions; 
        step::Int=1,
        oneletter=false,
        xlabel="Residue",
        ylabel="r / Å",
    )

Plot the contribution of each residue to the solute-solvent pair distribution function as a contour plot.
This function requires loading the `Plots` package.

# Arguments

- `rc::ResidueContributions`: The residue contributions to the solute-solvent pair distribution function,
   as computed by the `ResidueContributions` function.

# Optional arguments

- `step`: The step of the residue ticks in the x-axis of the plot. Default is 1 or will be set to show at most 20 ticks labels.
- `oneletter::Bool`: Use one-letter residue codes. Default is `false`. One-letter codes are only available for the 20 standard amino acids.
- `xlabel` and `ylabel`: Labels for the x and y axes. Default is `"Residue"` and `"r / Å"`.

This function will set color limits and choose the scale automatically. These parameters can be set
with the `clims` and `color` arguments of `Plots.contourf`. Other plot customizations can be 
done by passing other keyword arguments to this function, which will be passed to `Plots.contourf`.

# Example

```julia-repl
julia> using ComplexMixtures, Plots, PDBTools

julia> results = load("mddf.json")

julia> atoms = readPDB("system.pdb", "protein")

julia> rc = ResidueContributions(results, atoms; oneletter=true)

julia> plt = contourf(rc; step=5)
```

This will produce a plot with the contribution of each residue to the solute-solvent pair distribution function,
as a contour plot, with the residues in the x-axis and the distance in the y-axis.

To customize the plot, use the `Plot.contourf` keyword parameters, for example:

```julia-repl
julia> plt = contourf(rc; step=5, size=(800,400), title="Title", clims=(-0.1, 0.1))
```

!!! compat
    This function requires loading the `Plots` package and is available in
    ComplexMixtures v2.5.0 or greater.

    Support for all `Plots.contourf` parameters was introduced in ComplexMixtures v2.6.0.

"""
function Plots.contourf(
    rc::ResidueContributions;
    step::Union{Nothing,Integer}=nothing,
    oneletter=false,
    xlabel="Residue",
    ylabel="r / Å",
    clims=nothing,
    color=nothing,
    kargs...
)

    # Plot a contour curves with the density at each distance from each residue
    # colors, linewidths, etc. are defined here and can be tuned
    input_step = isnothing(step) ? 1 : step
    tick_range = firstindex(rc.xticks[1]):input_step:lastindex(rc.xticks[1])
    nticks = 50
    if isnothing(step) && length(tick_range) > nticks
        step = length(rc.resnums) ÷ nticks
        @warn """\n
            Consider using a step to set the number of residue ticks in the plot. 
            - step will be set to $step to display $nticks residue ticks.

        """ _line=nothing _file=nothing
        tick_range = first(tick_range):step:last(tick_range)
    end
    tick_marks = rc.xticks[1][tick_range]
    tick_labels = rc.xticks[2][tick_range]
    tick_resnums = rc.resnums[tick_range]
    warned = false
    xticks = if oneletter
        for i in eachindex(tick_labels)
            resnum = "$(tick_resnums[i])"
            # The following will allow custom residue names to be identified
            # in the tick label string, by removing the residue number
            residue_name = string(strip(tick_labels[i][1:(first(findlast(resnum, tick_labels[i]))-1)]))
            one_letter_code = PDBTools.oneletter(residue_name)
            if one_letter_code == "X"
                if !warned 
                    @warn "One-letter code for residue(s) not found. Using full residue name." _file=nothing _line=nothing
                    warned = true
                end
                one_letter_code = residue_name
            end
            tick_labels[i] = one_letter_code * resnum
        end
        (tick_marks, tick_labels)
    else
        (tick_marks, tick_labels)
    end
    clims, colorscale = _set_clims_and_colorscale!(rc; clims, colorscale=color)
    if colorscale == :tempo
        levels = 5
    else
        levels = 12
    end

    # density to plot
    rc_range = 1:length(rc)
    if length(rc_range) > 2000 
        rc_step = length(rc_range) ÷ 2000
        rc_range = 1:rc_step:length(rc)
        @warn """\n
            The number of residues to plot is too large. Will plot every $rc_step residues.

        """ _line=nothing _file=nothing
    end
    Plots.default(fontfamily="Computer Modern")
    plt = Plots.contourf(
        rc.xticks[1][rc_range],
        rc.d, hcat(rc[rc_range].residue_contributions...);
        color=Plots.cgrad(colorscale),
        linewidth=1,
        linecolor=:black,
        colorbar=:none,
        levels=levels,
        xlabel=xlabel,
        ylabel=ylabel,
        xticks=xticks, xrotation=60, xtickfont=Plots.font(8, "Computer Modern"),
        size=(700, 400),
        margin=0.5Plots.PlotMeasures.cm,
        framestyle=:box,
        clims=clims,
        kargs...
    )

    # return the plot
    return plt
end

#
# This is a legacy function, will probably be deprecated some day.
#
ComplexMixtures.contourf_per_residue(rc::ResidueContributions; kwargs...) = Plots.contourf(rc; kwargs...)
"""
    contourf_per_residue(
        results::Result, atoms::AbstractVector{<:PDBTools.Atom}; 
        residue_range=nothing, 
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
- `atoms::AbstractVector{<:Atom}`: The atoms of the solute.

# Optional arguments

- `residue_range`: The range of residues to plot. Default is all residues. Use
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
    results::Result, atoms::AbstractVector{<:Atom};
    residue_range::AbstractRange=resnum(first(atoms)):resnum(last(atoms)),
    dmin=1.5, dmax=3.5,
    type=:mddf,
    kargs...
)
    first_atom = findfirst(at -> resnum(at) == first(residue_range), atoms)
    last_atom = findfirst(at -> resnum(at) == last(residue_range), atoms)
    if isnothing(first_atom) || isnothing(last_atom)
        throw(ArgumentError("""\n

            The residue_range $residue_range is out of bounds.

            The first and last residue numbers of the selection are: $(resnum(first(atoms))) and $(resnum(last(atoms))).

        """
        ))
    end
    rc = ResidueContributions(results, atoms[first_atom:last_atom]; dmin, dmax, type)
    return Plots.contourf(rc; step=step(residue_range), kargs...)
end

@testitem "contourf/contourf_per_residue" begin
    using ComplexMixtures
    using ComplexMixtures.Testing
    using Plots
    using PDBTools
    # Load example output file (computed in the previous script)
    protein = readPDB(joinpath(Testing.data_dir, "Gromacs/system.pdb"), "protein")
    results = load(joinpath(Testing.data_dir, "Gromacs/protein_EMI.json"))
    tmpplot = tempname() * ".png"

    # Add some strange residue names to the protein
    for atom in protein
        if resname(atom) == "LEU"
            atom.resname = "lx"
        end
    end

    rc = ResidueContributions(results, select(protein, "residue >= 50 and residue <= 75"))
    plt = contourf(rc)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)

    # different color scale
    rc2 = copy(rc)
    rc2.residue_contributions[:, 10:end] .*= -1
    plt = contourf(rc2)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)
    rc3 = copy(rc)
    rc3.residue_contributions .*= -1
    plt = contourf(rc3)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)

    plt = contourf(rc)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)

    plt = contourf(rc; step=5)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)

    plt = contourf(rc; clims=(-1,1))
    savefig(plt, tmpplot)
    @test isfile(tmpplot)

    plt = contourf(rc; colorscale=:tempo)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)

    plt = contourf(rc; clims=(-1,1), colorscale=:tempo)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)

    # Slice and plot
    rc2 = rc[10:5:20]
    plt = contourf_per_residue(results, protein; residue_range=50:2:70)
    savefig(plt, tmpplot)
    @test isfile(tmpplot)

    plt = contourf_per_residue(results, protein; residue_range=50:75, oneletter=true)
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

end # module Plotting