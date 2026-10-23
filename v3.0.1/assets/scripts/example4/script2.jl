import Pkg;
Pkg.activate(".");

using PDBTools
using ComplexMixtures
using Plots
using LaTeXStrings

# Load a PDB file of the system
system = read_pdb("./equilibrated.pdb")

# Select the atoms corresponding to glycerol and water (using PDBTools)
glyc = select(system, "resname GLLM")
water = select(system, "water")

# Load previously computed mddfs
mddf_glyc = load("./mddf_glyc.json")
mddf_glyc_water = load("./mddf_glyc_water.json")

# Plot some group contributions to the MDDF. We select the atom names
# corresponding to each type of group of the glycerol molecule.  
hydroxyls = ["O1", "O2", "O3", "HO1", "HO2", "HO3"]
aliphatic = ["C1", "C2", "C3", "H11", "H12", "H2", "H31", "H32"]

#
# Extract the contributions of these atoms to the MDDFs
#
# glycerol-glycerol
mddf_glyc_hydroxyls = contributions(mddf_glyc, SoluteGroup(hydroxyls))
mddf_glyc_aliphatic = contributions(mddf_glyc, SoluteGroup(aliphatic))
# glycerol-water
mddf_glyc_water_hydroxyls = contributions(mddf_glyc_water, SoluteGroup(hydroxyls))
mddf_glyc_water_aliphatic = contributions(mddf_glyc_water, SoluteGroup(aliphatic))

#
# Plot the contributions
#
Plots.default(
    fontfamily="Computer Modern",
    linewidth=2,
    framestyle=:box,
    grid=false,
    label=nothing,
)
scalefontsizes();
scalefontsizes(1.2)
plot(layout=(2, 1))

#
# Group contributions to glycerol-glycerol auto correlation
#
plot!(mddf_glyc.d, mddf_glyc.mddf, label="Total", subplot=1)
plot!(mddf_glyc.d, mddf_glyc_hydroxyls, label="Hydroxyls", subplot=1)
plot!(mddf_glyc.d, mddf_glyc_aliphatic, label="Aliphatic", subplot=1)

#
# Group contributions to glycerol-water correlation
#
plot!(mddf_glyc_water.d, mddf_glyc_water.mddf, label="Total", subplot=2)
plot!(mddf_glyc_water.d, mddf_glyc_water_hydroxyls, label="Hydroxyls", subplot=2)
plot!(mddf_glyc_water.d, mddf_glyc_water_aliphatic, label="Aliphatic", subplot=2)

# plot settings
plot!(ylabel="Glyc-Glyc MDDF", xlim=(1.0, 8.0), subplot=1)
plot!(ylabel="Glyc-Water MDDF", xlim=(1.0, 8.0), subplot=2)
plot!(xlabel=L"\mathrm{Distance / \AA}", subplot=2)

# Save figure
savefig("./mddf_group_contributions.png")
println("Plot saved to mddf_group_contributions.png")
