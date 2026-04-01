import Pkg;
Pkg.activate(".");

using ComplexMixtures
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

# Load previusly saved results, computed in the previous script
results = load("./mddf.json")

# Plot with two subplots
plot(layout=(2, 1))

# Plot the total mddf
plot!(
    results.d,
    movavg(results.mddf, n=10).x, # Smooth example with a running average
    label="Total",
    subplot=1
)

# Plot DMF group contributions to the MDDF. We use a named tuple where
# the keys are the group names, and the values are the atom names of the group
groups = (
    "CO" => ["C", "O"], # carbonyl
    "N" => ["N"],
    "Methyl groups" => ["CC", "CT", "HC1", "HC2", "HC3", "HT1", "HT2", "HT3"],
)
for (group_label, group_atoms) in groups
    # Retrieve the contributions of the atoms of this group
    group_contrib = contributions(results, SolventGroup(group_atoms))
    # Plot the contributions of this groups, with the appropriate label
    plot!(
        results.d,
        movavg(group_contrib, n=10).x,
        label=group_label,
        subplot=1
    )
end

# Adjust scale and label of axis
plot!(xlim=(1, 5), ylabel="MDDF", subplot=1)

#
# Plot ACR group contributions to the MDDF. This is an interesting case,
# as the groups are repeated along the polymer chain
#
groups = (
    "CH_3" => ["CF", "HF1", "HF2", "HF3", "CL", "HL1", "HL2", "HL3"], # terminal methyles
    "CO" => ["OE1", "CD"], # carbonyl
    "NH_2" => ["NE2", "HE22", "HE21"], # amine
    "CHCH_2" => ["C", "H2", "H1", "CA", "HA"], # backbone
)
# Plot total mddf 
plot!(
    results.d,
    movavg(results.mddf, n=10).x, # Smooth example with a running average
    label="Total",
    subplot=2
)
# Plot group contributions
for (group_name, atom_names) in groups
    group_contrib = contributions(results, SoluteGroup(atom_names))
    plot!(
        results.d,
        movavg(group_contrib, n=10).x,
        label=latexstring("\\textrm{$group_name}"),
        subplot=2
    )
end

# Adjust scale and label of axis
plot!(
    xlim=(1, 5),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF", subplot=2
)
# Save figure
savefig("./mddf_groups.png")
println("Created figure file: ./mddf_groups.png")