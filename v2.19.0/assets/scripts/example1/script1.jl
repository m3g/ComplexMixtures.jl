# Activate environment in current directory
import Pkg;
Pkg.activate(".");

# Run this once, to install necessary packages:
# Pkg.add(["ComplexMixtures", "PDBTools", "Plots", "LaTeXStrings"])

# Load packages
using ComplexMixtures
using PDBTools
using Plots, Plots.Measures
using LaTeXStrings

# The complete trajectory file can be downloaded from (1Gb):
# https://www.dropbox.com/scl/fi/zfq4o21dkttobg2pqd41m/glyc50_traj.dcd?rlkey=el3k6t0fx6w5yiqktyx96gzg6&dl=0

# The example output file is available at:
# 
# Load PDB file of the system
atoms = read_pdb("./system.pdb")

# Select the protein and the GLYC molecules
protein = select(atoms, "protein")
glyc = select(atoms, "resname GLYC")

# Setup solute and solvent structures
solute = AtomSelection(protein, nmols=1)
solvent = AtomSelection(glyc, natomspermol=14)

# Path to the trajectory file
trajectory_file = "./glyc50_traj.dcd"

# Run mddf calculation, and save results
results = mddf(trajectory_file, solute, solvent, Options(bulk_range=(10.0, 15.0)))
save(results, "glyc50_results.json")
println("Results saved to glyc50_results.json")

#
# Produce plots
#
# Default options for plots 
Plots.default(
    fontfamily="Computer Modern",
    linewidth=2,
    framestyle=:box,
    label=nothing,
    grid=false
)

#
# The complete MDDF and the Kirkwood-Buff Integral
#
plot(layout=(1, 2))
# plot mddf
plot!(results.d, results.mddf,
    xlabel=L"r/\AA",
    ylabel="mddf",
    subplot=1
)
hline!([1], linestyle=:dash, linecolor=:gray, subplot=1)
# plot KB integral
plot!(results.d, results.kb / 1000, #to L/mol
    xlabel=L"r/\AA",
    ylabel=L"G_{us}/\mathrm{L~mol^{-1}}",
    subplot=2
)
# size and margin
plot!(size=(800, 300), margin=4mm)
savefig("./mddf.png")
println("Created plot mddf.png")

#
# Atomic contributions to the MDDF
#
hydroxyls = ["O1", "O2", "O3", "H1", "H2", "H3"]
aliphatic = ["C1", "C2", "HA", "HB", "HC", "HD"]
hydr_contrib = contributions(results, SolventGroup(hydroxyls))
aliph_contrib = contributions(results, SolventGroup(aliphatic))

plot(results.d, results.mddf,
    xlabel=L"r/\AA",
    ylabel="mddf",
    size=(600, 400)
)
plot!(results.d, hydr_contrib, label="Hydroxyls")
plot!(results.d, aliph_contrib, label="Aliphatic chain")
hline!([1], linestyle=:dash, linecolor=:gray)
savefig("./mddf_atom_contrib.png")
println("Created plot mddf_atom_contrib.png")