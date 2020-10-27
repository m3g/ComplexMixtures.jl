using Test
using ComplexMixtures, PDBTools
using Random
const CM = ComplexMixtures

@testset "NAMD" begin

  #
  # Tests with NAMD-DCD trajectory
  #
  
  dir="./data/NAMD"
  atoms = readPDB("$dir/structure.pdb")  
  
  # Example 1: protein-tmao
  
  R_save = CM.load("$dir/protein_tmao.json")
  protein = CM.Selection(select(atoms,"protein"),nmols=1)
  tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao) 
  options = CM.Options(stride=5,seed=1234567,nthreads=1)
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t
  
  # Example 2: water-tmao
  
  R_save = CM.load("$dir/water_tmao.json")
  tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  water = CM.Selection(select(atoms,"water"),natomspermol=3)
  traj = CM.Trajectory("$dir/trajectory.dcd",tmao,water) 
  options = CM.Options(stride=5,seed=1234567,nthreads=1)
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t
  
  # Example 3: tmao-tmao
  
  R_save = CM.load("$dir/tmao_tmao.json")
  tmao = CM.Selection(select(atoms,"resname TMAO"),natomspermol=14)
  traj = CM.Trajectory("$dir/trajectory.dcd",tmao) 
  options = CM.Options(stride=1,seed=1234567,nthreads=1)
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t
  
  # Example 3: water-water
  
  R_save = CM.load("$dir/water_water.json")
  water = CM.Selection(select(atoms,"water"),natomspermol=3)
  traj = CM.Trajectory("$dir/trajectory.dcd",water) 
  options = CM.Options(stride=5,seed=1234567,nthreads=1)
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t

end

@testset "Gromacs" begin

  #
  # Tests with NAMD-DCD trajectory
  #
  
  dir="./data/Gromacs"
  atoms = readPDB("$dir/system.pdb")  
  
  # Example 1: protein-EMIM
  
  R_save = CM.load("$dir/protein_EMI.json")
  protein = CM.Selection(select(atoms,"protein"),nmols=1)
  emi = CM.Selection(select(atoms,"resname EMI"),natomspermol=20)
  traj = CM.Trajectory("$dir/trajectory.xtc",protein,emi) 
  options = CM.Options(stride=5,seed=1234567,nthreads=1)
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t
  
  # Example 1: EMIM-DCA
  
  R_save = CM.load("$dir/EMI_DCA.json")
  emi = CM.Selection(select(atoms,"resname EMI"),natomspermol=20)
  dca = CM.Selection(select(atoms,"resname NC"),natomspermol=5)
  traj = CM.Trajectory("$dir/trajectory.xtc",emi,dca) 
  options = CM.Options(stride=5,seed=1234567,nthreads=1)
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t
  
  # Example 1: EMIM-EMIM
  
  R_save = CM.load("$dir/EMI_EMI.json")
  emi = CM.Selection(select(atoms,"resname EMI"),natomspermol=20)
  traj = CM.Trajectory("$dir/trajectory.xtc",emi) 
  options = CM.Options(stride=5,seed=1234567,nthreads=1)
  R = CM.mddf(traj,options)
  t = isapprox(R,R_save,debug=true) 
  @test t
  
end













