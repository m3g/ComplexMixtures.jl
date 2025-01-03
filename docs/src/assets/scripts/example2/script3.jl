import Pkg;
Pkg.activate(".");

using ComplexMixtures
using Plots
using EasyFit: movavg
using LaTeXStrings
using PDBTools

# Here we will produce a 2D plot of group contributions, splitting the
# contributions of each mer of the polymer into its chemical groups

# Chemical groups of the polymer monomers, defined by the atom types:
groups = (
    "CH_{3L}" => ["CF", "HF1", "HF2", "HF3"], # left-terminal methyl
    "CO" => ["OE1", "CD"], # carbonyl
    "NH_2" => ["NE2", "HE22", "HE21"], # amine 
    "CHCH_2" => ["C", "H2", "H1", "CA", "HA"], # backbone
    "CH_{3R}" => ["CL", "HL1", "HL2", "HL3"], # right-terminal methyles
)

system = readPDB("./equilibrated.pdb")
acr = select(system, "resname FACR or resname ACR or resname LACR")
results = load("./mddf.json")

# Here we split the polymer in residues, to extract the contribution of 
# each chemical group of each polymer mer independently
group_contribs = Vector{Float64}[]
labels = LaTeXString[]
for mer in eachresidue(acr)
    for (group_label, group_atoms) in groups
        if resnum(mer) != 1 && group_label == "CH_{3L}" # first residue only
            continue
        end
        if resnum(mer) != 5 && group_label == "CH_{3R}" # last residue only
            continue
        end
        # Filter the atoms of this mer that belong to the group
        mer_group_atoms = filter(at -> name(at) in group_atoms, mer)
        # Retrive the contribution of this mer atoms to the MDDF
        atoms_contrib = contributions(results, SoluteGroup(mer_group_atoms))
        # Smooth the contributions
        atoms_contrib = movavg(atoms_contrib; n=10).x
        # Add contributions to the group contributions list
        push!(group_contribs, atoms_contrib)
        # Push label to label list, using LaTeX formatting
        push!(labels, latexstring("\\textrm{$group_label}"))
    end
end

# Convert the group contributions to a matrix
group_contribs = stack(group_contribs)

# Find the indices of the limits of the map we want
idmin = findfirst(d -> d > 1.5, results.d)
idmax = findfirst(d -> d > 3.2, results.d)

# Plot contour map
Plots.default(fontfamily="Computer Modern")
contourf(
    1:length(labels),
    results.d[idmin:idmax],
    group_contribs[idmin:idmax, :],
    color=cgrad(:tempo), linewidth=1, linecolor=:black,
    colorbar=:none, levels=10,
    xlabel="Group", ylabel=L"r/\AA", xrotation=60,
    xticks=(1:length(labels), labels),
    margin=5Plots.Measures.mm # adjust margin 
)
savefig("./map2D_acr.png")
println("Plot saved to map2D_acr.png")