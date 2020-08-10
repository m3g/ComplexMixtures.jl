#
# Protein - Water
#

using PDBTools
using MDDF

atoms = PDBTools.readPDB("./structure.pdb")
solute_indexes = PDBTools.select(atoms,"protein")
solute = MDDF.Selection( solute_indexes, nmols=1 )
solvent_indexes = PDBTools.select(atoms,"resname TIP3")
solvent = MDDF.Selection( solvent_indexes, natomspermol=3 )
trajectory = MDDF.Trajectory("../NAMD/trajectory.dcd",solute,solvent)
options = MDDF.Options(binstep=0.2,n_random_samples=10)
@time R = MDDF.mddf_linkedcells_parallel(trajectory,options)

