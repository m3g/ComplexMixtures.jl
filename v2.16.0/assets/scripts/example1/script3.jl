import Pkg;
Pkg.activate(".");
using PDBTools
using ComplexMixtures

# PDB file of the system simulated
atoms = read_pdb("./system.pdb")

# Load results of a ComplexMixtures run
results = load("./glyc50_results.json")

# Compute the 3D density grid and output it to the PDB file
# here we use dmax=3.5 such that the the output file is not too large
grid = grid3D(results, atoms, "./grid.pdb"; dmin=1.5, dmax=3.5)