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
system = read_pdb("equilibrated.pdb")

# Load the pre-calculated MDDF of ethanol
mddf_ethanol_POPC = load("mddf_ethanol_POPC.json")

# Splitting the oleoyl chain into groups along the the chain. 
# The labels `CH_2` etc stand for `CH₂`, for example, in LaTeX notation, 
# for a nicer plot axis ticks formatting.
oleoyl_groups = (
    "CO" => ["O22", "C21"],
    "CH_2" => ["H2R", "H2S", "C22"],
    "CH_2" => ["C23", "H3R", "H3S"],
    "CH_2" => ["C24", "H4R", "H4S"],
    "CH_2" => ["C25", "H5R", "H5S"],
    "CH_2" => ["C26", "H6R", "H6S"],
    "CH_2" => ["C27", "H7R", "H7S"],
    "CH_2" => ["C28", "H8R", "H8S"],
    "CH" => ["C29", "H91"],
    "CH" => ["C210", "H101"],
    "CH_2" => ["C211", "H11R", "H11S"],
    "CH_2" => ["C212", "H12R", "H12S"],
    "CH_2" => ["C213", "H13R", "H13S"],
    "CH_2" => ["C214", "H14R", "H14S"],
    "CH_2" => ["C215", "H15R", "H15S"],
    "CH_2" => ["C216", "H16R", "H16S"],
    "CH_2" => ["C217", "H17R", "H17S"],
    "CH_3" => ["C218", "H18R", "H18S", "H18T"]
)

# Format tick labels with LaTeX
labels_o = [latexstring("\\textrm{$key}") for (key, val) in oleoyl_groups]

# We first collect the contributions of each group into a vector of vectors:
gcontrib = Vector{Float64}[] # empty vector of vectors
for (group_name, group_atoms) in oleoyl_groups
    group_contributions = contributions(mddf_ethanol_POPC, SoluteGroup(group_atoms))
    push!(gcontrib, movavg(group_contributions; n=10).x)
end

# Convert the vector of vectors into a matrix
gcontrib = stack(gcontrib)

# Find the indices of the MDDF where the distances are between 1.5 and 3.0 Å
idmin = findfirst(d -> d > 1.5, mddf_ethanol_POPC.d)
idmax = findfirst(d -> d > 3.0, mddf_ethanol_POPC.d)

# The plot will have two lines, the first plot will contain the 
# oleoyl groups contributions, and the second plot will contain the
# contributions of the palmitoyl groups.
plot(layout=(2, 1))

# Plot the contributions of the oleoyl groups
contourf!(
    1:length(oleoyl_groups),
    mddf_ethanol_POPC.d[idmin:idmax],
    gcontrib[idmin:idmax, :],
    color=cgrad(:tempo), linewidth=1, linecolor=:black,
    colorbar=:none, levels=10,
    ylabel=L"r/\AA", xrotation=60,
    xticks=(1:length(oleoyl_groups), labels_o),
    subplot=1,
)
annotate!(14, 2.7, text("Oleoyl", :left, 12, plot_font), subplot=1)

#
# Repeat procedure for the palmitoyl groups
#
palmitoyl_groups = (
    "CO" => ["C31", "O32"],
    "CH_2" => ["C32", "H2X", "H2Y"],
    "CH_2" => ["C33", "H3X", "H3Y"],
    "CH_2" => ["C34", "H4X", "H4Y"],
    "CH_2" => ["C35", "H5X", "H5Y"],
    "CH_2" => ["C36", "H6X", "H6Y"],
    "CH_2" => ["C37", "H7X", "H7Y"],
    "CH_2" => ["C38", "H8X", "H8Y"],
    "CH_2" => ["C39", "H9X", "H9Y"],
    "CH_2" => ["C310", "H10X", "H10Y"],
    "CH_2" => ["C311", "H11X", "H11Y"],
    "CH_2" => ["C312", "H12X", "H12Y"],
    "CH_2" => ["C313", "H13X", "H13Y"],
    "CH_2" => ["C314", "H14X", "H14Y"],
    "CH_2" => ["C315", "H15X", "H15Y"],
    "CH_3" => ["C316", "H16X", "H16Y", "H16Z"],
)

# Format tick labels with LaTeX
labels_p = [latexstring("\\textrm{$key}") for (key, val) in palmitoyl_groups]

# We first collect the contributions of each group into a # vector of vectors:
gcontrib = Vector{Float64}[] # empty vector of vectors
for (group_name, group_atoms) in palmitoyl_groups
    group_contributions = contributions(mddf_ethanol_POPC, SoluteGroup(group_atoms))
    push!(gcontrib, movavg(group_contributions; n=10).x)
end

# Convert the vector of vectors into a matrix
gcontrib = stack(gcontrib)

# Plot the contributions of the palmitoyl groups
contourf!(
    1:length(palmitoyl_groups),
    mddf_ethanol_POPC.d[idmin:idmax],
    gcontrib[idmin:idmax, :],
    color=cgrad(:tempo), linewidth=1, linecolor=:black,
    colorbar=:none, levels=10,
    xlabel="Group",
    ylabel=L"r/\AA", xrotation=60,
    xticks=(1:length(labels_p), labels_p),
    bottom_margin=0.5Plots.Measures.cm,
    subplot=2,
)
annotate!(12, 2.7, text("Palmitoyl", :left, 12, plot_font), subplot=2)

savefig("POPC_ethanol_chains.png")
println("The plot was saved as POPC_ethanol_chains.png")

#
# Now, plot a similar map for the water interactions with the POPC chain
#
mddf_water_POPC = load("mddf_water_POPC.json")
plot(layout=(2, 1))

# Plot the contributions of the oleoyl groups
# We first collect the contributions of each group into a vector of vectors:
gcontrib = Vector{Float64}[] # empty vector of vectors
for (group_name, group_atoms) in oleoyl_groups
    group_contributions = contributions(mddf_water_POPC, SoluteGroup(group_atoms))
    push!(gcontrib, movavg(group_contributions; n=10).x)
end
# Convert the vector of vectors into a matrix
gcontrib = stack(gcontrib)
# Plot matrix as density map
contourf!(
    1:length(oleoyl_groups),
    mddf_water_POPC.d[idmin:idmax],
    gcontrib[idmin:idmax, :],
    color=cgrad(:tempo), linewidth=1, linecolor=:black,
    colorbar=:none, levels=10,
    ylabel=L"r/\AA", xrotation=60,
    xticks=(1:length(oleoyl_groups), labels_o), subplot=1,
)
annotate!(14, 2.7, text("Oleoyl", :left, 12, plot_font), subplot=1)

# Plot the contributions of the palmitoyl groups
# We first collect the contributions of each group into a vector of vectors:
gcontrib = Vector{Float64}[] # empty vector of vectors
for (group_name, group_atoms) in palmitoyl_groups
    group_contributions = contributions(mddf_water_POPC, SoluteGroup(group_atoms))
    push!(gcontrib, movavg(group_contributions; n=10).x)
end
# Convert the vector of vectors into a matrix
gcontrib = stack(gcontrib)
# Plot matrix as density map
contourf!(
    1:length(palmitoyl_groups),
    mddf_water_POPC.d[idmin:idmax],
    gcontrib[idmin:idmax, :],
    color=cgrad(:tempo), linewidth=1, linecolor=:black,
    colorbar=:none, levels=10,
    xlabel="Group",
    ylabel=L"r/\AA", xrotation=60,
    xticks=(1:length(palmitoyl_groups), labels_o),
    subplot=2,
    bottom_margin=0.5Plots.Measures.cm,
)
annotate!(12, 2.7, text("Palmitoyl", :left, 12, plot_font), subplot=2)

savefig("POPC_water_chains.png")
println("The plot was saved as POPC_water_chains.png")