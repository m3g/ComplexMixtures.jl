# Activate environment in current directory
import Pkg;
Pkg.activate(".");

# Run this once, to install necessary packages:
# Pkg.add(["ComplexMixtures", "PDBTools", "Plots", "LaTeXStrings"])

# Load packages
using ComplexMixtures
using PDBTools
using Plots

# Load PDB file of the system
atoms = readPDB("./system.pdb")

# Select the protein and the GLYC molecules
protein = select(atoms, "protein")

# Load example output file (computed in the previous script)
example_output = "./glyc50_results.json"
results = load(example_output)

#
# Plot a 2D map showing the contributions of some residues
# - the residue range referes to serial residue indices in the
#   structure provided.
#
contourf_per_residue(
  results, protein;
  residue_range=70:110,
  dmin=1.5, dmax=3.5,
  oneletter=true
)
savefig("./density2D.png")
println("Created plot density2D.png")
