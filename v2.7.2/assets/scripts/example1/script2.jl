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
#
residue_contributions = ResidueContributions(
  results, 
  select(protein, "resnum >= 70 and resnum <= 110")
)
contourf(residue_contributions; oneletter=true)
savefig("./density2D.png")
println("Created plot density2D.png")
