import Pkg;
Pkg.activate(".");

using PDBTools
using ComplexMixtures
using Plots
using LaTeXStrings
using EasyFit: movavg

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
x = mddf_glyc.d # distances
# Total mddf
y = movavg(mddf_glyc.mddf; n=10).x
plot!(x, y, label="Total", subplot=1)
# Hydroxyls
y = movavg(mddf_glyc_hydroxyls; n=10).x
plot!(x, y, label="Hydroxyls", subplot=1)
# Aliphatic
y = movavg(mddf_glyc_aliphatic; n=10).x
plot!(x, y, label="Aliphatic", subplot=1)

#
# Group contributions to glycerol-water correlation
#
x = mddf_glyc_water.d # distances
# Total mddf
y = movavg(mddf_glyc_water.mddf; n=10).x
plot!(x, y, label="Total", subplot=2)
# Hydroxyls
y = movavg(mddf_glyc_water_hydroxyls; n=10).x
plot!(x, y, label="Hydroxyls", subplot=2)
# Aliphatic
y = movavg(mddf_glyc_water_aliphatic; n=10).x
plot!(x, y, label="Aliphatic", subplot=2)

# plot settings
plot!(ylabel="Glyc-Glyc MDDF", xlim=(1.0, 8.0), subplot=1)
plot!(ylabel="Glyc-Water MDDF", xlim=(1.0, 8.0), subplot=2)
plot!(xlabel=L"\mathrm{Distance / \AA}", subplot=2)

# Save figure
savefig("./mddf_group_contributions.png")
println("Plot saved to mddf_group_contributions.png")
