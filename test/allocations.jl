
@testitem "Allocations" begin
  using BenchmarkTools
  using ComplexMixtures.Testing
  using ComplexMixtures, PDBTools
  const CM = ComplexMixtures

  dir = "$(Testing.data_dir)/NAMD"
  atoms = readPDB("$dir/structure.pdb")

  options = Options(lastframe=1, nthreads=1, silent=true, seed=321, StableRNG=true)
  t_options = @benchmark Options(lastframe=1, seed=321, StableRNG=true, nthreads=1, silent=true) samples = 1 evals = 1
  @test t_options.allocs == 0

  protein = Selection(select(atoms, "protein"), nmols=1)
  t_selection1 = @benchmark Selection(select($atoms, "protein"), nmols=1) samples = 1 evals = 1
  @test t_selection1.allocs == 51

  tmao = Selection(select(atoms, "resname TMAO"), natomspermol=14)
  t_selection2 = @benchmark Selection(select($atoms, "resname TMAO"), natomspermol=14) samples = 1 evals = 1
  @test t_selection2.allocs == 186139

  trajfile = "$dir/trajectory.dcd" # because of the interpolation of @benchmark
  traj = Trajectory(trajfile, protein, tmao)
  t_trajectory = @benchmark Trajectory($trajfile, $protein, $tmao) samples = 1 evals = 1
  @test t_trajectory.allocs == 844

  R = Result(traj, options)
  t_result = @benchmark Result($traj, $options) samples = 1 evals = 1
  @test t_result.allocs == 87

  CM.opentraj!(traj)
  CM.firstframe!(traj)
  t_nextframe = @benchmark CM.nextframe!($traj) samples = 1 evals = 1
  @test t_nextframe.allocs == 35

  RNG = CM.init_random(options)
  t_RNG = @benchmark CM.init_random($options) samples = 1 evals = 1
  @test t_RNG.allocs == 2

  system = CM.setup_PeriodicSystem(traj, options)
  buff = CM.Buffer(traj, R)
  @. buff.solute_read = traj.x_solute
  @. buff.solvent_read = traj.x_solvent
  CM.update_unitcell!(system, CM.getsides(traj, 1))
  t_mddf_frame = @benchmark CM.mddf_frame!($R, $system, $buff, $options, $RNG) samples = 1 evals = 1
  @test t_mddf_frame.allocs < 100

end
