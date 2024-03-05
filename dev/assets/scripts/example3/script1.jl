import Pkg;
Pkg.activate(".");

using PDBTools
using ComplexMixtures
using Plots
using LaTeXStrings
using EasyFit: movavg

# The full trajectory file is available at: 
# https://www.dropbox.com/scl/fi/hcenxrdf8g8hfbllyakhy/traj_POPC.dcd?rlkey=h9zivtwgya3ivva1i6q6xmr2p&dl=0
trajectory_file = "./traj_POPC.dcd"

# Load a PDB file of the system
system = readPDB("./equilibrated.pdb")

# Select the atoms corresponding to glycerol and water
popc = select(system, "resname POPC")
water = select(system, "water")
ethanol = select(system, "resname ETOH")

# Set the complete membrane as the solute. We use nmols=1 here such
# that the membrane is considered a single solute in the calculation. 
solute = AtomSelection(popc, nmols=1)

# Compute water-POPC distribution and KB integral 
solvent = AtomSelection(water, natomspermol=3)

# Set the trajectory structure
trajectory = Trajectory(trajectory_file, solute, solvent)

# We want to get reasonably converged KB integrals, which usually
# require large solute domains. Distribution functions converge 
# rapidly (~10Angs or less), on the other side.
options = Options(bulk_range=(20.0, 25.0))

# Compute the mddf and associated properties
mddf_water_POPC = mddf(trajectory, options)

# Save results to file for later use
save(mddf_water_POPC, "./mddf_water_POPC.json")
println("Results saved to ./mddf_water_POPC.json file")

# Compute ethanol-POPC distribution and KB integral 
solvent = AtomSelection(ethanol, natomspermol=9)
traj = Trajectory(trajectory_file, solute, solvent)
mddf_ethanol_POPC = mddf(traj, options)

# Save results for later use
save(mddf_ethanol_POPC, "./mddf_ethanol_POPC.json")
println("Results saved to ./mddf_ethanol_POPC.json file")

#
# Plot the MDDF and KB integrals
#
# Plot defaults
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=false,
    palette=:tab10
)
scalefontsizes(); scalefontsizes(1.3)

#
# Plots cossolvent-POPC MDDFs in subplot 1
#
plot(layout=(2,1))
# Water MDDF
plot!(
    mddf_water_POPC.d, # distances
    movavg(mddf_water_POPC.mddf,n=10).x, # water MDDF - smoothed
    label="Water",
    subplot=1
)
# Ethanol MDDF
plot!(
    mddf_ethanol_POPC.d, # distances
    movavg(mddf_ethanol_POPC.mddf,n=10).x, # water MDDF - smoothed
    label="Ethanol",
    subplot=1
)
# Plot settings
plot!(
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF",
    xlim=(0,10),
    subplot=1
)

#
# Plot cossolvent-POPC KB integrals in subplot 2
#
# Water KB
plot!(
    mddf_water_POPC.d, # distances
    mddf_water_POPC.kb, # water KB
    label="Water",
    subplot=2
)
# Ethanol KB
plot!(
    mddf_ethanol_POPC.d, # distances
    mddf_ethanol_POPC.kb, # ethanol KB
    label="Ethanol",
    subplot=2
)
# Plot settings
plot!(
    xlabel=L"\textrm{Distance / \AA}",
    ylabel=L"\textrm{KB~/~L~mol^{-1}}",
    xlim=(0,10),
    subplot=2
)

savefig("popc_water_ethanol_mddf_kb.png")
println("Plot saved to popc_water_ethanol_mddf_kb.png file")


