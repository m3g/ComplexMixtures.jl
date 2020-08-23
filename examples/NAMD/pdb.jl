#
# Protein - TMAO (compare new and old implementations)
#

include("../../src/MDDF.jl")
using PDBTools

atoms = PDBTools.readPDB("./structure.pdb")

protein = PDBTools.select(atoms,"protein")
solute = MDDF.Selection(protein,nmols=1 )

water = PDBTools.select(atoms,"water")
solvent = MDDF.Selection(water,natomspermol=3)

options = MDDF.Options()
trajectory = MDDF.PDBTraj("./traj200.pdb",solute,solvent)

pdb = MDDF.mddf(trajectory,options)

