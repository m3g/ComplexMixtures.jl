module ComplexMixturesPlotting

using ComplexMixtures
using PDBTools: Residue, residue_ticks
import Plots

function Plots.contourf(
    results::Result, residues::AbstractVector{Residue};
    residue_range=1:length(residues),
    dmin = 1.5, dmax = 3.5,
)

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
    xticks = residue_ticks(protein, first=first(residue_range), last=last(residue_range))
    
    # Plot a contour courves with the density at each distance from
    # each residue
    Plots.default(fontfamily="Computer Modern")
    plt = contourf(residue_range, results.d[idmin:idmax], rescontrib[idmin:idmax, residue_range],
      color=cgrad(:tempo), linewidth=1, linecolor=:black,
      colorbar=:none, levels=5,
      xlabel="Residue", ylabel="r / Å",
      xticks=xticks, xrotation=60,
      xtickfont=font(8, "Computer Modern"),
      size=(700, 400),
      margin=0.5Plots.PlotMeasures.cm
    )

    return plt
end

end