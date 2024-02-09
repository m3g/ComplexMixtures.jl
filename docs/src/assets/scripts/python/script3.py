# Load packages
import ComplexMixtures as cm
import matplotlib.pyplot as plt

# Read the pdb file and set solvent and solute groups
atoms = cm.readPDB("./system.pdb")
protein = cm.select(atoms, "protein")
glyc = cm.select(atoms, "resname GLYC")

# load previusly computed MDDF results
results = cm.load("./glyc50.json")

# Select atoms by name
hydroxyls = cm.list(["O1","O2","O3","H1","H2","H3"])
aliphatic = cm.list(["C1","C2","HA","HB","HC","HD"])

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