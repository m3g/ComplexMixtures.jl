import Pkg;
Pkg.activate(".");

using ComplexMixtures
using PDBTools
using Plots
using LaTeXStrings
using EasyFit: movavg

# Some default settings for the plots
plot_font = "Computer Modern"
Plots.default(
    fontfamily=plot_font,
    linewidth=1.5,
    framestyle=:box,
    label=nothing,
    grid=false,
)
scalefontsizes();
scalefontsizes(1.3);

# Read system PDB file
system = readPDB("equilibrated.pdb")
ethanol = select(system, "resname ETOH")

# Load the pre-calculated MDDF of ethanol
mddf_ethanol_POPC = load("mddf_ethanol_POPC.json")

#
# Contributions of the ethanol groups
#
# Define the groups using selections. Set a dict, in which the keys are the group names
# and the values are the selections
groups = (
    "Hydroxyl" => select(ethanol, "name O1 or name HO1"),
    "Aliphatic chain" => select(ethanol, "not name O1 and not name HO1"),
)
# plot the total mddf and the contributions of the groups
x = mddf_ethanol_POPC.d
plot(x, movavg(mddf_ethanol_POPC.mddf, n=10).x, label="Total MDDF")
for (group_name, group_atoms) in groups
    cont = contributions(mddf_ethanol_POPC, SolventGroup(group_atoms))
    y = movavg(cont, n=10).x
    plot!(x, y, label=group_name)
end
# Plot settings
plot!(
    xlim=(1, 8),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF"
)
savefig("./mddf_ethanol_groups.png")
println("The plot was saved as mddf_ethanol_groups.png")