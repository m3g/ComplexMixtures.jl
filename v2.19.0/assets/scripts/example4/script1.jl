import Pkg;
Pkg.activate(".");

using PDBTools
using ComplexMixtures
using Plots
using LaTeXStrings
using EasyFit: movavg

# Load a PDB file of the system
system = read_pdb("./equilibrated.pdb")

# The full trajectory file is available at: 
# https://www.dropbox.com/scl/fi/ag7k2d7i9d7ivbd5zmtl9/traj_Glyc.dcd?rlkey=93i31a5ytlzb34ulzjz315eyq&dl=0
trajectory_file = "./traj_Glyc.dcd"

# Select the atoms corresponding to glycerol and water (using PDBTools)
glyc = select(system, "resname GLLM")
water = select(system, "water")

# Compute Glycerol-Glycerol auto correlation mddf: first we initialize the 
# AtomSelection object, informing the number of atoms per molecule of Glycerol
glyc_selection = AtomSelection(glyc, natomspermol=14)

# We define a large solute domain (large dbulk) to obtain a good convergence
# for the KB integral. The mddf converges at much shorter distances.   
options = Options(bulk_range=(20.0, 25.0))
mddf_glyc = mddf(trajectory_file, glyc_selection, options)

# Save results for later analysis
save(mddf_glyc, "./mddf_glyc.json")

# Compute water-glycerol mddf: glyc_selection is the solute, and water is the solvent
solvent = AtomSelection(water, natomspermol=3)
mddf_glyc_water = mddf(trajectory_file, glyc_selection, solvent, options)

# Save results for later analysis
save(mddf_glyc_water, "./mddf_glyc_water.json")

#
# Plot the MDDFs 
#
Plots.default(
    fontfamily="Computer Modern",
    linewidth=2,
    framestyle=:box,
    grid=false,
    label=nothing,
)
scalefontsizes();
scalefontsizes(1.3)
plot(layout=(2, 1))

# glycerol-glycerol auto correlation
x = mddf_glyc.d # distances
y = movavg(mddf_glyc.mddf, n=10).x # the mddf (using movavg to smooth noise)
plot!(x, y, label="Glycerol-Glycerol", subplot=1)

# water-glycerol correlation
x = mddf_glyc_water.d
y = movavg(mddf_glyc_water.mddf, n=10).x
plot!(x, y, label="Glycerol-Water", subplot=1)
plot!(ylabel="MDDF", xlim=(1.5, 8), subplot=1)

# Plot the KB integrals
# glycrerol-glycerol
y = movavg(mddf_glyc.kb, n=10).x
plot!(x, y, subplot=2)

# water-glycerol
y = movavg(mddf_glyc_water.kb, n=10).x
plot!(x, y, subplot=2)

# plot settings
plot!(
    xlabel=L"\textrm{Distance / \AA}",
    ylabel=L"\textrm{KB~/~cm^3~mol^{-1}}",
    xlim=(0, 20),
    subplot=2
)

# Save plot
savefig("./mddf_kb.png")
println("Plot saved to mddf_kb.png")

