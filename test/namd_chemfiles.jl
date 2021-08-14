using Test
using ComplexMixtures, PDBTools
const CM = ComplexMixtures

@testset "NAMD with ChemFiles" begin

  #
  # Tests with NAMD-DCD trajectory
  #
  
  dir="./data/NAMD"
  atoms = readPDB("$dir/structure.pdb")  
  options = Options(stride=5,seed=321,StableRNG=true,nthreads=1,silent=true)
  
  # Example 1: protein-tmao
  
  # save(R,"$dir/protein_tmao.json")
  R_save = load("$dir/protein_tmao.json")
  protein = Selection(select(atoms,"protein"),nmols=1)
  tmao = Selection(select(atoms,"resname TMAO"),natomspermol=14)
  traj = CM.Trajectory("$dir/trajectory.dcd",protein,tmao,chemfiles=true) 
  R = mddf(traj,options)
  @test isapprox(R,R_save,debug=true) 
  
  # Test save and load
  save(R,"test.json")
  R_load = load("test.json")
  @test R_load ≈ R_save
  #rm("./test.json",force=true)

  # Example 2: water-tmao
  
  # save(R,"$dir/water_tmao.json")
  R_save = load("$dir/water_tmao.json")
  tmao = Selection(select(atoms,"resname TMAO"),natomspermol=14)
  water = Selection(select(atoms,"water"),natomspermol=3)
  traj = CM.Trajectory("$dir/trajectory.dcd",tmao,water,chemfiles=true) 
  R = mddf(traj,options)
  @test isapprox(R,R_save,debug=true) 
  
  # Example 3: tmao-tmao
  
  # save(R,"$dir/tmao_tmao.json")
  R_save = load("$dir/tmao_tmao.json")
  tmao = Selection(select(atoms,"resname TMAO"),natomspermol=14)
  traj = CM.Trajectory("$dir/trajectory.dcd",tmao,tmao,chemfiles=true) 
  R = mddf(traj,options)
  @test isapprox(R,R_save,debug=true) 
  
  # Example 3: water-water
  
  # save(R,"$dir/water_water.json")
  R_save = load("$dir/water_water.json")
  water = Selection(select(atoms,"water"),natomspermol=3)
  traj = CM.Trajectory("$dir/trajectory.dcd",water,water,chemfiles=true) 
  R = mddf(traj,options)
  @test isapprox(R,R_save,debug=true) 

end
