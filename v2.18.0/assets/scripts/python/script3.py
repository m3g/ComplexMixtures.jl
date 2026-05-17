# Load packages
import ComplexMixtures as cm
import matplotlib.pyplot as plt

# Read the pdb file and set solvent and solute groups
atoms = cm.read_pdb("./system.pdb")
protein = cm.select(atoms, "protein")
glyc = cm.select(atoms, "resname GLYC")

# load previously computed MDDF results
results = cm.load("./glyc50.json")

# Select GLYC subgroups by atoms by name
hydroxyls = cm.select(glyc, "name O1 O2 O3 H1 H2 H3")
aliphatic = cm.select(glyc, "name C1 C2 HA HB HC HD")

# Extract the contributions of the Glycerol hydroxyls and aliphatic groups
hydr_contributions = cm.contributions(results, cm.SolventGroup(hydroxyls))
aliph_contributions = cm.contributions(results, cm.SolventGroup(aliphatic))

# Plot
plt.plot(results.d, results.mddf, label="Total MDDF")
plt.plot(results.d, hydr_contributions, label="Hydroxyls")
plt.plot(results.d, aliph_contributions, label="Aliphatic")
plt.legend()
plt.xlabel("distance / Angs")
plt.ylabel("MDDF")
plt.savefig("group_contributions.png")

#
# Extract the contributions of specific residues in the protein
#
my_selection = cm.select(atoms, "protein and residue 2 5 7")
residue_contributions = cm.contributions(results, cm.SoluteGroup(my_selection))
plt.plot(results.d, results.mddf, label="Total MDDF")
plt.plot(results.d, hydr_contributions, label="Residue contributions")
plt.legend()
plt.xlabel("distance / Angs")
plt.ylabel("MDDF")
plt.savefig("residue_contributions.png")

#
# Note: to select more complex groups of atoms, you can use the `select_with_vmd` function,
# which allows you to use VMD-style selection strings. For using so, you need to have
# VMD installed and available in your PATH, or provide it with the `vmd` keyword.
# 
# See: https://m3g.github.io/PDBTools.jl/stable/selections/#Using-VMD
#