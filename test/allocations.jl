using Test
using ComplexMixtures, PDBTools
const CM = ComplexMixtures

@testset "Allocations" begin

  if Sys.islinux() && VERSION == v"1.5.3" 

    dir="./data/NAMD"
    atoms = readPDB("$dir/structure.pdb")  

    options = CM.Options(lastframe=1,nthreads=1,silent=true,seed=321,StableRNG=true)
    t_options = @allocated CM.Options(lastframe=1,seed=321,StableRNG=true,nthreads=1,silent=true)
    @test t_options == 0

    protein = CM.Selection(select(atoms,"protein"),nmols=1)
    t_selection1 = @allocated CM.Selection(select(atoms,"protein"),nmols=1)
    @test t_selection1 == 3064544

    tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
    t_selection2 = @allocated CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
    @test t_selection2 == 7080880

    traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
    t_trajectory = @allocated CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
    @test abs(t_trajectory - 660240) == 240

    samples = CM.Samples(md=(traj.solvent.nmols-1)/2,random=options.n_random_samples)
    t_samples = @allocated CM.Samples(md=(traj.solvent.nmols-1)/2,random=options.n_random_samples)
    @test t_samples == 256

    R = CM.Result(traj,options)
    t_result = @allocated CM.Result(traj,options)
    @test t_result == 5968736

    framedata = CM.FrameData(traj,R)
    t_framedata = @allocated CM.FrameData(traj,R)
    @test t_framedata == 228528

    CM.nextframe!(traj)
    t_nextframe = @allocated CM.nextframe!(traj)
    @test t_nextframe == 624

    RNG = CM.init_random(options)
    t_RNG = @allocated CM.init_random(options) 
    @test t_RNG == 0

    CM.mddf_frame!(1,framedata,options,RNG,R)
    t_mddf_frame = @allocated CM.mddf_frame!(1,framedata,options,RNG,R)
    @test t_mddf_frame == 160

  end
  @test true

end

