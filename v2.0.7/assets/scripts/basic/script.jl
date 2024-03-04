# Activate environment (see the Installation -> Recommended Workflow manual section)
import Pkg;
Pkg.activate(".");

# Load packages
using ComplexMixtures
using PDBTools
using Plots

# Load PDB file of the system
atoms = readPDB("./system.pdb")

# Select the protein and the TMAO molecules
protein = select(atoms, "protein")
tmao = select(atoms, "resname TMAO")

# Setup solute and solvent structures. We need to provide
# either the number of atoms per molecule, or the number
# of molecules in each selection.
solute = AtomSelection(protein, nmols=1)
solvent = AtomSelection(tmao, natomspermol=14)

# Setup the Trajectory structure: this will define which
# coordinates are used to compute the MDDF when reading
# the trajectory file.
trajectory = Trajectory("./trajectory.dcd", solute, solvent)

# Run the calculation and get results: this is the computationally
# intensive part of the calculation.
results = mddf(trajectory)

# Save the results to recover them later if required
save(results, "./results.json")

# Plot the some of the most important results.
#
# - The results.d array contains the distances. 
# - The results.mddf array contains the MDDF.
# - The results.kb array contains the Kirkwood-Buff integrals.
#
plot(results.d, results.mddf, xlabel="d / Å", ylabel="MDDF") # plot the MDDF
savefig("./mddf.pdf")
plot(results.d, results.kb, xlabel="d / Å", ylabel="KB / cm³ mol⁻¹") # plot the KB 
savefig("./kb.pdf")