using Revise
using ComplexMixtures, PDBTools
const CM = ComplexMixtures

dir="./test/data/NAMD"
atoms = readPDB("$dir/structure.pdb")  

@time options = Options(lastframe=1,seed=321,StableRNG=true,nthreads=1,silent=true)

@time tmao = Selection(select(atoms,"resname TMAO"),natomspermol=14)

# cross-correlation

@time protein = Selection(select(atoms,"protein"),nmols=1)

@time traj = Trajectory("$dir/trajectory.dcd",protein,tmao) 

@time samples = CM.Samples(md=(traj.solvent.nmols-1)/2,random=options.n_random_samples)

@time R = Result(traj,options)

@time framedata = CM.FrameData(traj,R)

@time CM.nextframe!(traj)

@time RNG = CM.init_random(options) 

@time CM.mddf_frame!(1,framedata,options,RNG,R)

# self

@time traj = Trajectory("$dir/trajectory.dcd",tmao) 

@time samples = CM.Samples(md=(traj.solvent.nmols-1)/2,random=options.n_random_samples)

@time R = Result(traj,options)

@time framedata = CM.FrameData(traj,R)

@time CM.nextframe!(traj)

@time RNG = CM.init_random(options) 

@time CM.mddf_frame_self!(1,framedata,options,RNG,R)



