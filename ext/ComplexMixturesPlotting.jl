module ComplexMixturesPlotting

import ComplexMixtures
using ComplexMixtures: Result, SoluteGroup, contributions
using PDBTools: Residue, residue_ticks, Atom, eachresidue
import Plots

#
# Function to plot the contour map fo the contribution of each residue
# to the solute-solvent pair distribution function
#
function ComplexMixtures.contourf_per_residue(
    results::Result, atoms::AbstractVector{Atom};
    residue_range=1:length(eachresidue(atoms)),
    dmin=1.5, dmax=3.5,
    oneletter=false,
)

    # collect the list of residues (using PDBTools)
    residues = collect(eachresidue(atoms))

    # of bins of the mddf histogram (length(results.d)) and a number of 
    # columns equal to the number of residues
    rescontrib = zeros(length(results.d), length(residues))

    # Each column is then filled up with the contributions of each residue
    for (ires, residue) in enumerate(residues)
        rescontrib[:, ires] .= contributions(results, SoluteGroup(residue))
    end

    # Plot only for distances within 1.5 and 3.5:
    idmin = findfirst(d -> d > dmin, results.d)
    idmax = findfirst(d -> d > dmax, results.d)

    # Obtain pretty labels for the residues in the x-axis
    xticks = residue_ticks(atoms, first=first(residue_range), last=last(residue_range); oneletter)

    # Plot a contour courves with the density at each distance from
    # each residue
    Plots.default(fontfamily="Computer Modern")
    plt = Plots.contourf(residue_range, results.d[idmin:idmax], rescontrib[idmin:idmax, residue_range],
        color=Plots.cgrad(:tempo), linewidth=1, linecolor=:black,
        colorbar=:none, levels=5,
        xlabel="Residue", ylabel="r / Å",
        xticks=xticks, xrotation=60,
        xtickfont=Plots.font(8, "Computer Modern"),
        size=(700, 400),
        margin=0.5Plots.PlotMeasures.cm
    )

    return plt
end

end