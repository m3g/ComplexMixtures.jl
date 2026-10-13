# The ComplexMixtures.py file is assumed to be in the current
# directory.
# Obtain it from: 
# https://m3g.github.io/ComplexMixtures.jl/stable/assets/ComplexMixtures.py 
import ComplexMixtures as cm

# Load the pdb file of the system using `PDBTools`:
atoms = cm.read_pdb("./system.pdb")

# Create arrays of atoms with the protein and Glycerol atoms, 
# using the `select` function of the `PDBTools` package:
#
# Note: for more complex selections, you can use the `cm.select_with_vmd` function,
# which allows you to use VMD-style selection strings. For using so, you need to have
# VMD installed and available in your PATH, or provide it with the `vmd` keyword.
# See: https://m3g.github.io/PDBTools.jl/stable/selections/#Using-VMD
#
protein = cm.select(atoms,"protein")
glyc = cm.select(atoms,"resname GLYC")
water = cm.select(atoms,"water")

# Setup solute and solvent structures, required for computing the MDDF, 
# with `AtomSelection` function of the `ComplexMixtures` package:
solute = cm.AtomSelection(protein, nmols=1)
solvent = cm.AtomSelection(glyc, natomspermol=14)

# Run the calculation and get results:
results = cm.mddf("./glyc50_sample.dcd", solute, solvent, cm.Options(bulk_range=(10.0,12.0)))

# Save the reults to recover them later if required
cm.save(results,"./glyc50.json")
print("Results saved to glyc50.json")

# Compute the water distribution function around the protein:
solvent = cm.AtomSelection(water, natomspermol=3)
results = cm.mddf("./glyc50_sample.dcd", solute, solvent, cm.Options(bulk_range=(10.0,12.0)))
cm.save(results,"./water.json")
print("Results saved to water.json")