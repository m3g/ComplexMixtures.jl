import Pkg;
Pkg.activate(".");

using PDBTools
using ComplexMixtures
using Plots
using LaTeXStrings

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

# We define a large solute domain (large cutoff) to obtain a good convergence
# for the KB integral. The mddf converges at much shorter distances.   
options = Options(cutoff=25.0)
mddf_glyc = mddf(trajectory_file, glyc_selection, options)

# Save results for later analysis
save("./mddf_glyc.json", mddf_glyc)

# Compute water-glycerol mddf: glyc_selection is the solute, and water is the solvent
solvent = AtomSelection(water, natomspermol=3)
mddf_glyc_water = mddf(trajectory_file, glyc_selection, solvent, options)

# Save results for later analysis
save("./mddf_glyc_water.json", mddf_glyc_water)

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
plot!(mddf_glyc.d, mddf_glyc.mddf, label="Glycerol-Glycerol", subplot=1)

# water-glycerol correlation
plot!(mddf_glyc_water.d, mddf_glyc_water.mddf, label="Glycerol-Water", subplot=1)
plot!(ylabel="MDDF", xlim=(1.5, 8), subplot=1)

# plot KB integrals
plot!(
    mddf_glyc_water.d, 
    [ mddf_glyc.kb mddf_glyc_water.kb ],
    xlabel=L"\textrm{Distance / \AA}",
    ylabel=L"\textrm{KB~/~cm^3~mol^{-1}}",
    xlim=(0, 20),
    subplot=2
)

# Save plot
savefig("./mddf_kb.png")
println("Plot saved to mddf_kb.png")

