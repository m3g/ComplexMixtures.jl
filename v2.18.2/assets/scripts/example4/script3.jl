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

# 2D maps plot of group contributions

# Glycerol groups
groups = (
    "OH" => ["O1", "HO1"], # first hydroxyl
    "CH_2" => ["C1", "H11", "H12"], # first CH2
    "OH" => ["O2", "HO2"], # second hydroxyl
    "CH" => ["C2", "H2"], # CH
    "CH_2" => ["C3", "H31", "H32"], # second CH2
    "OH" => ["O3", "HO3"] # third hydroxyl
)
labels = latexstring.("\\textrm{$name}" for (name, atoms) in groups)

#
# Contributions to Glycerol-Glycerol autocorrelation
# First, create a vector of vectors, in which each component carries the
# contributions of each Glycerol group to the MDDF
#
group_contrib = Vector{Float64}[] # empty vector of vectors
for (name, atoms) in groups
    push!(group_contrib, contributions(mddf_glyc, SolventGroup(atoms)))
end

# Convert output to a matrix to plot a 2D map
group_contrib = stack(group_contrib)

# The distance range to plot
idmin = findfirst(d -> d > 1.5, mddf_glyc.d)
idmax = findfirst(d -> d > 3.0, mddf_glyc.d)

#
# plot the map
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

# map of contributions to the Glycerol-Glycerol autocorrelation
contourf!(
    1:length(groups),
    mddf_glyc.d[idmin:idmax],
    group_contrib[idmin:idmax, :],
    color=cgrad(:tempo), linewidth=1, linecolor=:black,
    colorbar=:none, levels=5,
    xticks=(1:length(groups), labels), xrotation=60,
    ylabel=L"r/\AA",
    subplot=1
)

# Water-glycerol interactions (Glycerol contributions)
group_contrib = Vector{Float64}[] # empty vector of vectors
for (name, atoms) in groups
    push!(group_contrib, contributions(mddf_glyc_water, SoluteGroup(atoms)))
end

# Convert output to a matrix to plot a 2D map
group_contrib = stack(group_contrib)

# map of the contributions of Glycerol groups to the Glycerol-Water correlation
contourf!(
    1:length(groups),
    mddf_glyc_water.d[idmin:idmax],
    group_contrib[idmin:idmax, :],
    color=cgrad(:tempo), linewidth=1, linecolor=:black,
    colorbar=:none, levels=5,
    xticks=(1:length(groups), labels), xrotation=60,
    ylabel=L"r/\AA",
    subplot=2
)
plot!(
    xlabel="Glycerol group",
    bottommargin=0.5Plots.Measures.cm,
    subplot=2
)

savefig("./GlycerolWater_map.png")
println("Plot saved to GlycerolWater_map.png")