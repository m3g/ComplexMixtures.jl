# Activate environment in current directory
import Pkg; Pkg.activate(".")

# Run this once, to install necessary packages:
# Pkg.add(["ComplexMixtures", "PDBTools", "Plots", "LaTeXStrings"])

# Load packages
using ComplexMixtures
using PDBTools
using Plots, Plots.Measures
using LaTeXStrings

# The complete trajectory file can be downloaded from (3Gb):
# https://drive.google.com/file/d/14M30jDHRwUM77hzbDphgbu8mcWFBcQrX/view?usp=sharing

# The example output file is available at:
# 
# Load PDB file of the system
atoms = readPDB("./system.pdb")

# Select the protein and the GLYC molecules
protein = select(atoms, "protein")
glyc = select(atoms, "resname GLYC")

# Load example output file (computed in the previous script)
example_output = "./glyc50_results.json"
results = load(example_output)

#
# Plot a 2D map showing the contributions of some residues
#
residues = collect(eachresidue(protein))

# We will plot only the range 70:110, for clarity
irange = 70:110

# We create matrix of with a number of rows equal to the number
# of bins of the mddf histogram (length(results.d)) and a number of 
# columns equal to the number of residues
rescontrib = zeros(length(results.d), length(residues))

# Each column is then filled up with the contributions of each residue
for (ires, residue) in enumerate(residues)
    rescontrib[:, ires] .= contributions(results, SoluteGroup(residue))
end

# Plot only for distances within 1.5 and 3.5:
idmin = findfirst(d -> d > 1.5, results.d)
idmax = findfirst(d -> d > 3.5, results.d)

# Obtain pretty labels for the residues in the x-axis
xticks = PDBTools.residue_ticks(protein, first=70, last=110)

# Plot a contour courves with the density at each distance from
# each residue
Plots.default(fontfamily="Computer Modern")
contourf(irange, results.d[idmin:idmax], rescontrib[idmin:idmax, irange],
  color=cgrad(:tempo), linewidth=1, linecolor=:black,
  colorbar=:none, levels=5,
  xlabel="Residue", ylabel=L"r/\AA",
  xticks=xticks, xrotation=60,
  xtickfont=font(8, "Computer Modern"),
  size=(700, 400),
  margin=0.5Plots.PlotMeasures.cm
)
savefig("./density2D.png")
println("Created plot density2D.png")