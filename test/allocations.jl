using Test
using ComplexMixtures, PDBTools
using Random
const CM = ComplexMixtures

@testset "Allocations" begin

  dir="./data/NAMD"
  atoms = readPDB("$dir/structure.pdb")  

  options = CM.Options(lastframe=1,seed=1234567,nthreads=1,silent=true)
  t_options = @allocated options = CM.Options(lastframe=1,seed=1234567,nthreads=1,silent=true)
  @test t_options == 0

  protein = CM.Selection(select(atoms,"protein"),nmols=1)
  t_selection1 = @allocated protein = CM.Selection(select(atoms,"protein"),nmols=1)
  @test t_selection1 == 3064544

  tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  t_selection2 = @allocated tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  @test t_selection2 == 7080880

  traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
  t_trajectory = @allocated traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
  @test abs(t_trajectory - 660240) < 124

  samples = CM.Samples(md=(traj.solvent.nmols-1)/2,random=options.n_random_samples)
  t_samples = @allocated samples = CM.Samples(md=(traj.solvent.nmols-1)/2,random=options.n_random_samples)
  @test t_samples == 256

  R = CM.Result(traj,options)
  t_result = @allocated R = CM.Result(traj,options)
  @test t_result == 5968688

  framedata = CM.FrameData(traj,R)
  t_framedata = @allocated framedata = CM.FrameData(traj,R)
  @test t_framedata == 228528

  CM.nextframe!(traj)
  t_nextframe = @allocated CM.nextframe!(traj)
  @test t_nextframe == 624

  iframe = 1
  CM.mddf_frame!(iframe,framedata,options,R)
  t_mddf_frame = @allocated CM.mddf_frame!(iframe,framedata,options,R)
  @test t_mddf_frame == 144

end

