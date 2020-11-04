using Revise
using ComplexMixtures, PDBTools
using Random
const CM = ComplexMixtures

dir="./test/data/NAMD"
atoms = readPDB("$dir/structure.pdb")  

@time options = CM.Options(lastframe=1,seed=1234567,nthreads=1,silent=true)

@time tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)

# cross-correlation

@time protein = CM.Selection(select(atoms,"protein"),nmols=1)

@time traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 

@time samples = CM.Samples(md=(traj.solvent.nmols-1)/2,random=options.n_random_samples)

@time R = CM.Result(traj,options)

@time framedata = CM.FrameData(traj,R)

@time CM.nextframe!(traj)

iframe = 1
@time CM.mddf_frame!(iframe,framedata,options,R)

# self

@time traj = CM.Trajectory("$dir/trajectory.dcd",tmao) 

@time samples = CM.Samples(md=(traj.solvent.nmols-1)/2,random=options.n_random_samples)

@time R = CM.Result(traj,options)

@time framedata = CM.FrameData(traj,R)

@time CM.nextframe!(traj)

iframe = 1
@time CM.mddf_frame_self!(iframe,framedata,options,R)



