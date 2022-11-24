import Pkg; Pkg.activate(".")

using ComplexMixtures
using PDBTools

data_dir = "$(@__DIR__)/../data/NAMD"

atoms = readPDB(data_dir*"/structure.pdb")
options = Options(stride=1,seed=321,StableRNG=true,nthreads=1,silent=true)
protein = Selection(select(atoms, "protein"), nmols=1)
tmao = Selection(select(atoms, "resname TMAO"), natomspermol=14)
traj = Trajectory("$(data_dir)/trajectory.dcd", protein, tmao)
R1 = mddf(traj, options)

#water = Selection(select(atoms, "water"), natomspermol=3)
#traj = Trajectory("$(data_dir)/trajectory.dcd", protein, water)
#R2 = mddf(traj, options)
#
#atoms = readPDB(data_dir*"/../simple.pdb")
#options = Options(stride=1,seed=321,StableRNG=true,nthreads=1,silent=true)
#protein = Selection(select(atoms, "protein"), nmols=1)
#water = Selection(select(atoms, "resname WAT"), natomspermol=3)
#traj = Trajectory("$(data_dir)/../simple.pdb", protein, water, format="PDBTraj")
#R3 = mddf(traj, options)
