using Test
using BenchmarkTools
using ComplexMixtures, PDBTools
const CM = ComplexMixtures

@testset "Allocations" begin

  dir="./data/NAMD"
  atoms = readPDB("$dir/structure.pdb")  

  options = Options(lastframe=1,nthreads=1,silent=true,seed=321,StableRNG=true)
  t_options = @benchmark Options(lastframe=1,seed=321,StableRNG=true,nthreads=1,silent=true) samples=1 evals=1
  @test t_options.allocs == 0

  protein = Selection(select(atoms,"protein"),nmols=1)
  t_selection1 = @benchmark Selection(select($atoms,"protein"),nmols=1) samples=1 evals=1
  @test t_selection1.allocs == 62082

  tmao = Selection(select(atoms,"resname TMAO"),natomspermol=14)
  t_selection2 = @benchmark Selection(select($atoms,"resname TMAO"),natomspermol=14) samples=1 evals=1
  @test t_selection2.allocs == 186144

  trajfile = "$dir/trajectory.dcd" # because of the interpolation of @benchmark
  traj = Trajectory(trajfile,protein,tmao) 
  t_trajectory = @benchmark Trajectory($trajfile,$protein,$tmao) samples=1 evals=1
  @test t_trajectory.allocs == 839

  samples = CM.Samples(md=(traj.solvent.nmols-1)/2,random=options.n_random_samples)
  t_samples = @benchmark CM.Samples(md=($traj.solvent.nmols-1)/2,random=$options.n_random_samples) samples=1 evals=1
  @test t_samples.allocs == 0

  R = Result(traj,options)
  t_result = @benchmark Result($traj,$options) samples=1 evals=1
  @test t_result.allocs == 66

  framedata = CM.FrameData(traj,R)
  t_framedata = @benchmark CM.FrameData($traj,$R) samples=1 evals=1
  @test t_framedata.allocs == 206

  CM.nextframe!(traj)
  t_nextframe = @benchmark CM.nextframe!($traj) samples=1 evals=1
  @test t_nextframe.allocs == 34

  RNG = CM.init_random(options)
  t_RNG = @benchmark CM.init_random($options) samples=1 evals=1
  @test t_RNG.allocs == 2

  CM.mddf_frame!(1,framedata,options,RNG,R)
  t_mddf_frame = @benchmark CM.mddf_frame!(1,$framedata,$options,$RNG,$R) samples=1 evals=1
  @test t_mddf_frame.allocs == 0

end

