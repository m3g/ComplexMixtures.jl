
include("../src/MDDF.jl")

# Here we use the PDBTools package to read the pdb file (from http://github.com/m3g/PDBTools)
using PDBTools
atoms = PDBTools.readPDB("./structure.pdb")

# The solute is water, with 3 atoms
solute_indexes = [ atom.index for atom in filter( atom -> atom.resname == "TIP3", atoms ) ]
solute = MDDF.Solute( solute_indexes, natomspermol=3 )

# The solvent is TMAO, which has 14 atoms. Use the natomspermol to indicate how many
# atoms each molecule has, such that there is no ambiguity on how to split the coordinates 
# of the selection into individual molecules.
solvent_indexes = [ atom.index for atom in filter( atom -> atom.resname == "TMAO", atoms ) ]
solvent = MDDF.Solvent( solvent_indexes, natomspermol=14 )

# Initialize trajectroy data structure and open input stream
trajectory = MDDF.NamdDCD("./trajectory.dcd",solute,solvent)

# Input options for the calcualtion
options = MDDF.Options(output="example.dat")

# Run MDDF calculation, and get the resutls in the R structure
R = MDDF.mddf_naive(trajectory,options)

include("./plots.jl")

