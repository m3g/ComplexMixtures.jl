using Test
using ComplexMixtures, PDBTools
using Random
const CM = ComplexMixtures

@testset "Internals" begin

  dir="./data/NAMD"
  atoms = readPDB("$dir/structure.pdb")  

  options = CM.Options(lastframe=1,seed=1234567,nthreads=1,silent=true)
  t_options = @timed options = CM.Options(lastframe=1,seed=1234567,nthreads=1,silent=true)
  @test t_options.bytes == 0

  protein = CM.Selection(select(atoms,"protein"),nmols=1)
  t_selection1 = @timed protein = CM.Selection(select(atoms,"protein"),nmols=1)
  @test t_selection1.bytes == 3064544

  tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  t_selection2 = @timed tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  @test t_selection2.bytes == 7080880

  traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
  t_trajectory = @timed traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
  @test t_trajectory.bytes == 660240

  samples = CM.Samples(md=(traj.solvent.nmols-1)/2,random=options.n_random_samples)
  t_samples = @timed samples = CM.Samples(md=(traj.solvent.nmols-1)/2,random=options.n_random_samples)
  @test t_samples.bytes == 256

  R = CM.Result(traj,options)
  t_result = @timed R = CM.Result(traj,options)
  @test t_result.bytes == 5968352

  framedata = CM.FrameData(traj,R)
  t_framedata = @timed framedata = CM.FrameData(traj,R)
  @test t_framedata.bytes == 229984

  CM.nextframe!(traj)
  t_nextframe = @timed CM.nextframe!(traj)
  @test t_nextframe.bytes == 624

  iframe = 1
  CM.mddf_frame!(iframe,framedata,options,R)
  t_mddf_frame = @timed CM.mddf_frame!(iframe,framedata,options,R)
  @test t_mddf_frame.bytes == 144

end

@testset "NAMD" begin

  #
  # Tests with NAMD-DCD trajectory
  #
  
  dir="./data/NAMD"
  atoms = readPDB("$dir/structure.pdb")  
  options = CM.Options(stride=5,seed=1234567,nthreads=1,silent=true)
  
  # Example 1: protein-tmao
  
  # CM.save(R,"$dir/protein_tmao.json")
  R_save = CM.load("$dir/protein_tmao.json")
  protein = CM.Selection(select(atoms,"protein"),nmols=1)
  tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t
  
  # Test save and load
  CM.save(R,"test.json")
  R_load = CM.load("test.json")
  @test R_load ≈ R_save

  # Example 2: water-tmao
  
  # CM.save(R,"$dir/water_tmao.json")
  R_save = CM.load("$dir/water_tmao.json")
  tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  water = CM.Selection(select(atoms,"water"),natomspermol=3)
  traj = CM.Trajectory("$dir/trajectory.dcd",tmao,water) 
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t
  
  # Example 3: tmao-tmao
  
  # CM.save(R,"$dir/tmao_tmao.json")
  R_save = CM.load("$dir/tmao_tmao.json")
  tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  traj = CM.Trajectory("$dir/trajectory.dcd",tmao) 
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t
  
  # Example 3: water-water
  
  # CM.save(R,"$dir/water_water.json")
  R_save = CM.load("$dir/water_water.json")
  water = CM.Selection(select(atoms,"water"),natomspermol=3)
  traj = CM.Trajectory("$dir/trajectory.dcd",water) 
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t

end

@testset "Gromacs" begin

  #
  # Tests with Gromacs-XSC trajectory
  #
  
  dir="./data/Gromacs"
  atoms = readPDB("$dir/system.pdb")  
  options = CM.Options(stride=5,seed=1234567,nthreads=1,silent=true)
  
  # Example 1: protein-EMIM
  
  # CM.save(R,"$dir/protein_EMI.json")
  R_save = CM.load("$dir/protein_EMI.json")
  protein = CM.Selection(select(atoms,"protein"),nmols=1)
  emi = CM.Selection(select(atoms,"resname EMI"),natomspermol=20)
  traj = CM.Trajectory("$dir/trajectory.xtc",protein,emi) 
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t
  
  # Example 1: EMIM-DCA
  
  # CM.save(R,"$dir/EMI_DCA.json")
  R_save = CM.load("$dir/EMI_DCA.json")
  emi = CM.Selection(select(atoms,"resname EMI"),natomspermol=20)
  dca = CM.Selection(select(atoms,"resname NC"),natomspermol=5)
  traj = CM.Trajectory("$dir/trajectory.xtc",emi,dca) 
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t
  
  # Example 1: EMIM-EMIM
  
  # CM.save(R,"$dir/EMI_EMI.json")
  R_save = CM.load("$dir/EMI_EMI.json")
  emi = CM.Selection(select(atoms,"resname EMI"),natomspermol=20)
  traj = CM.Trajectory("$dir/trajectory.xtc",emi) 
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t
  
end


@testset "PDB" begin

  #
  # Tests with trajectory given in a PDB file
  #
  
  dir="./data/PDB"
  atoms = readPDB("$dir/trajectory.pdb","model 1")  
  options = CM.Options(stride=1,seed=1234567,nthreads=1,silent=true)
 
  # Example 1: protein-tmao
  
  # CM.save(R,"$dir/protein_tmao.json")
  R_save = CM.load("$dir/protein_tmao.json")
  protein = CM.Selection(select(atoms,"protein"),nmols=1)
  tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  traj = CM.Trajectory("$dir/trajectory.pdb",protein,tmao,format="PDBTraj") 
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t
  
  # Example 3: tmao-tmao
  
  # CM.save(R,"$dir/tmao_tmao.json")
  R_save = CM.load("$dir/tmao_tmao.json")
  tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  traj = CM.Trajectory("$dir/trajectory.pdb",tmao,format="PDBTraj") 
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t
  
end 

@testset "Merge" begin

  dir="./data/NAMD"
  atoms = readPDB("$dir/structure.pdb")  
  tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  water = CM.Selection(select(atoms,"water"),natomspermol=3)
  
  # CM.save(R,"$dir/merged.json")
  R_save = CM.load("$dir/merged.json")

  options = CM.Options(firstframe=1,lastframe=2,seed=1234567,nthreads=1,silent=true)
  traj = CM.Trajectory("$dir/trajectory.dcd",tmao,water) 
  R1 = CM.mddf(traj,options)

  options = CM.Options(firstframe=3,lastframe=6,seed=1234567,nthreads=1,silent=true)
  traj = CM.Trajectory("$dir/trajectory.dcd",tmao,water) 
  R2 = CM.mddf(traj,options)

  R = CM.merge([R1,R2])
  t = isapprox(R,R_save,debug=true) 
  @test t

end
