#
# Protein - Water (no plots, to compute external timming)
#

using MDDF
using PDBTools

# Here we use the PDBTools package to read the pdb file (from http://github.com/m3g/PDBTools)
atoms = PDBTools.readPDB("./structure.pdb")

# The solute is a single protein molecule (infinte dilution case). In this case,
# use the option nmols=1
solute_indexes = PDBTools.select(atoms,"protein")
solute = MDDF.Selection( solute_indexes, nmols=1 )

# The solvent is Water, which has 3 atoms. Use the natomspermol to indicate how many
# atoms each molecule has, such that there is no ambiguity on how to split the coordinates 
# of the selection into individual molecules.
solvent_indexes = PDBTools.select(atoms,"water")
solvent = MDDF.Selection( solvent_indexes, natomspermol=3 )

# Initialize trajectroy data structure and open input stream
trajectory = MDDF.Trajectory("./trajectory.dcd",solute,solvent)

# Input options for the calcualtion
options = MDDF.Options(binstep=0.2)

# Run MDDF calculation, and get the resutls in the R structure
@time R = MDDF.mddf_linkedcells_parallel(trajectory,options)

