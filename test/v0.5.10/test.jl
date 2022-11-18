import Pkg; Pkg.activate(".")

using ComplexMixtures
using PDBTools

data_dir = "$(@__DIR__)/../data/NAMD"

atoms = readPDB(data_dir*"/structure.pdb")
options = Options(stride=5,seed=321,StableRNG=true,nthreads=1,silent=true)
protein = Selection(select(atoms, "protein"), nmols=1)
tmao = Selection(select(atoms, "resname TMAO"), natomspermol=14)
traj = Trajectory("$(data_dir)/trajectory.dcd", protein, tmao)
R = mddf(traj, options)


