import Pkg;
Pkg.activate(".");

using PDBTools
using ComplexMixtures
using Plots
using LaTeXStrings
using EasyFit: movavg

# The full trajectory file is available at: 
# https://www.dropbox.com/scl/fi/jwafhgxaxuzsybw2y8txd/traj_Polyacry.dcd?rlkey=p4bn65m0pkuebpfm0hf158cdm&dl=0 
trajectory_file = "./traj_Polyacry.dcd"

# Load a PDB file of the system
system = read_pdb("./equilibrated.pdb")

# Select the atoms corresponding DMF molecules
dmf = select(system, "resname DMF")

# Select the atoms corresponding to the Poly-acrylamide
acr = select(system, "resname FACR or resname ACR or resname LACR")

# Set the solute and the solvent selections for ComplexMixtures
solute = AtomSelection(acr, nmols=1)
solvent = AtomSelection(dmf, natomspermol=12)

# Use a large dbulk distance for better KB convergence
options = Options(bulk_range=(20.0, 25.0))

# Compute the mddf and associated properties
results = mddf(trajectory_file, solute, solvent, options)

# Save results to file for later use
save(results, "./mddf.json")
println("Results saved to ./mddf.json file")

# Plot the MDDF and KB integrals
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=1.5,
    framestyle=:box,
    label=nothing,
    grid=false,
    palette=:tab10
)
scalefontsizes();
scalefontsizes(1.3);

# Plot the MDDF of DMF relative to PolyACR and its corresponding KB integral
plot(layout=(2, 1))
plot!(
    results.d,
    movavg(results.mddf, n=9).x, # Smooth example with a running average
    ylabel="MDDF",
    xlims=(0, 20),
    subplot=1,
)

# Plot the KB integral
plot!(
    results.d,
    movavg(results.kb, n=9).x, # smooth kb
    xlabel=L"\textrm{Distance / \AA}",
    ylabel=L"\textrm{KB~/~cm^2~mol^{-1}}",
    xlim=(-1, 20),
    subplot=2
)
savefig("./mddf_kb.png")
println("Plot saved to mddf_kb.png")