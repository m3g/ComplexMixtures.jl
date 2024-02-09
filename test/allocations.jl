
@testitem "Allocations" begin
    using BenchmarkTools
    using ComplexMixtures.Testing
    using ComplexMixtures, PDBTools

    dir = "$(Testing.data_dir)/NAMD"
    atoms = readPDB("$dir/structure.pdb")

    options = Options(
        lastframe = 1,
        nthreads = 1,
        silent = true,
        seed = 321,
        StableRNG = true,
    )
    t = @benchmark Options(
        lastframe = 1,
        seed = 321,
        StableRNG = true,
        nthreads = 1,
        silent = true,
    ) samples = 1 evals = 1
    @test t.allocs == 0

    protein = AtomSelection(select(atoms, "protein"), nmols = 1)
    t_selection1 =
        @benchmark AtomSelection(select($atoms, "protein"), nmols = 1) samples = 1 evals = 1
    @test t_selection1.allocs < 100 

    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol = 14)
    t_selection2 =
        @benchmark AtomSelection(select($atoms, "resname TMAO"), natomspermol = 14) samples = 1 evals = 1
    @test t_selection2.allocs < 200000

    trajfile = "$dir/trajectory.dcd" # because of the interpolation of @benchmark
    traj = Trajectory(trajfile, protein, tmao)
    t_trajectory =
        @benchmark Trajectory($trajfile, $protein, $tmao) samples = 1 evals = 1
    @test t_trajectory.allocs < 1000

    R = Result(traj, options)
    t_result = @benchmark Result($traj, $options) samples = 1 evals = 1
    @test t_result.allocs < 5000

    ComplexMixtures.opentraj!(traj)
    ComplexMixtures.firstframe!(traj)
    t_nextframe = @benchmark ComplexMixtures.nextframe!($traj) samples = 1 evals = 1
    @test t_nextframe.allocs < 100

    RNG = ComplexMixtures.init_random(options)
    t_RNG = @benchmark ComplexMixtures.init_random($options) samples = 1 evals = 1
    @test t_RNG.allocs < 5

    system = ComplexMixtures.setup_PeriodicSystem(traj, options)
    buff = ComplexMixtures.Buffer(traj, R)
    @. buff.solute_read = traj.x_solute
    @. buff.solvent_read = traj.x_solvent
    ComplexMixtures.update_unitcell!(system, ComplexMixtures.convert_unitcell(ComplexMixtures.getunitcell(traj)))
    t_mddf_frame =
        @benchmark ComplexMixtures.mddf_frame!($R, $system, $buff, $options, 1.0, $RNG) samples = 1 evals = 1
    @test t_mddf_frame.allocs < 100

end
