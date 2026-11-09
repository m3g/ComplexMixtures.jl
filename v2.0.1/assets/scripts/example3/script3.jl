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
    linewidth=2,
    framestyle=:box,
    label=nothing,
    grid=false,
)
scalefontsizes();
scalefontsizes(1.3);

# Read system PDB file
system = readPDB("equilibrated.pdb")

# Load the pre-calculated MDDF of water
mddf_water_POPC = load("mddf_water_POPC.json")

#
# Here we define the POPC groups, from the atom names. Each group
# is a vector of atom names, and the keys are the group names.
#
groups = Dict(
    "Choline" => ["N", "C12", "H12A", "H12B", "C13", "H13A", "H13B", "H13C", "C14",
        "H14A", "H14B", "H14C", "C15", "H15A", "H15B", "H15C", "C11", "H11A", "H11B"],
    "Phosphate" => ["P", "O13", "O14", "O12"],
    "Glycerol" => ["O11", "C1", "HA", "HB", "C2", "HS", "O21", "C3", "HX", "HY", "O31"],
    "Oleoyl" => ["O22", "C21", "H2R", "H2S", "C22", "C23", "H3R", "H3S", "C24", "H4R", "H4S",
        "C25", "H5R", "H5S", "C26", "H6R", "H6S", "C27", "H7R", "H7S", "C28", "H8R", "H8S",
        "C29", "H91", "C210", "H101", "C211", "H11R", "H11S", "C212", "H12R", "H12S",
        "C213", "H13R", "H13S", "C214", "H14R", "H14S", "C215", "H15R", "H15S",
        "C216", "H16R", "H16S", "C217", "H17R", "H17S", "C218", "H18R", "H18S", "H18T"],
    "Palmitoyl" => ["C31", "O32", "C32", "H2X", "H2Y", "C33", "H3X", "H3Y", "C34", "H4X", "H4Y",
        "C35", "H5X", "H5Y", "C36", "H6X", "H6Y", "C37", "H7X", "H7Y", "C38", "H8X",
        "H8Y", "C39", "H9X", "H9Y", "C310", "H10X", "H10Y", "C311", "H11X", "H11Y",
        "C312", "H12X", "H12Y", "C313", "H13X", "H13Y", "C314", "H14X", "H14Y", "C315",
        "H15X", "H15Y", "C316", "H16X", "H16Y", "H16Z"],
)

#
# plot the total mddf and the contributions of the groups
#
x = mddf_water_POPC.d
plot(
    x,
    movavg(mddf_water_POPC.mddf, n=10).x,
    label="Total water-POPC MDDF"
)
for (group_name, group_atoms) in pairs(groups)
    cont = contributions(mddf_water_POPC, SoluteGroup(group_atoms))
    y = movavg(cont, n=10).x
    plot!(x, y, label=group_name)
end
# Plot settings
plot!(
    xlim=(1, 5),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF"
)
savefig("./mddf_POPC_water_groups.png")
println("The plot was saved as mddf_POPC_water_groups.png")

