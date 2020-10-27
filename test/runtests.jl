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
  # Tests with NAMD-DCD trajectory
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













