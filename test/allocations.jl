
@testitem "Allocations" begin
    using PDBTools
    using ComplexMixtures
    using ComplexMixtures.Testing
    using BenchmarkTools

    @kwdef struct Allocs
        prodbuild::Bool = haskey(ENV, "BUILD_IS_PRODUCTION_BUILD") && ENV["BUILD_IS_PRODUCTION_BUILD"] == "true"
        allocs::Int
    end
    Allocs(allocs::Int) = Allocs(; allocs)
    import Base: ==, >, <
    ==(a::Int, b::Allocs) = b.prodbuild ? a == b.allocs : true
    <(a::Int, b::Allocs) = b.prodbuild ? a < b.allocs : true
    ==(a::Allocs, b::Int) = a.prodbuild ? a.allocs == b : true
    <(a::Allocs, b::Int) = a.prodbuild ? a.allocs < b : true


    dir = "$(Testing.data_dir)/NAMD"
    atoms = readPDB("$dir/structure.pdb")

    options = Options(
        lastframe=1,
        nthreads=1,
        silent=true,
        seed=321,
        StableRNG=true,
    )
    t = @benchmark Options(
        lastframe=1,
        seed=321,
        StableRNG=true,
        nthreads=1,
        silent=true,
    ) samples = 1 evals = 1
    @test t.allocs == Allocs(0) 

    prot_atoms = select(atoms, "protein")
    indices = index.(prot_atoms)
    b = @benchmark AtomSelection($(Int.(indices)), nmols=1, natomspermol=11, group_names=$(String[]), group_atom_indices=$(Vector{Int}[])) samples = 1 evals = 1
    @test b.allocs == Allocs(0)

    protein = AtomSelection(prot_atoms, nmols=1)
    t_selection1A = @benchmark AtomSelection($prot_atoms, nmols=1) samples = 1 evals = 1
    @test t_selection1A.allocs < Allocs(1500) # one String per atom  name
    n = String.(name.(prot_atoms))
    t_selection1B = @benchmark AtomSelection($prot_atoms, nmols=1, group_names=$n, group_atom_indices=$(Vector{Int}[])) samples = 1 evals = 1
    @test t_selection1B.allocs <= Allocs(5)

    tmao_atoms = select(atoms, "resname TMAO")
    tmao = AtomSelection(tmao_atoms, natomspermol=14)
    t_selection2 =
        @benchmark AtomSelection($tmao_atoms, natomspermol=14) samples = 1 evals = 1
    @test t_selection2.allocs <= Allocs(100)

    trajfile = "$dir/trajectory.dcd" # because of the interpolation of @benchmark
    traj = Trajectory(trajfile, protein, tmao)
    t_trajectory =
        @benchmark Trajectory($trajfile, $protein, $tmao) samples = 1 evals = 1
    @test t_trajectory.allocs < Allocs(1000)

    R = Result(traj, options)
    t_result = @benchmark Result($traj, $options) samples = 1 evals = 1
    @test t_result.allocs < Allocs(5000)

    ComplexMixtures.opentraj!(traj)
    ComplexMixtures.firstframe!(traj)
    t_nextframe = @benchmark ComplexMixtures.nextframe!($traj) samples = 1 evals = 1
    @test t_nextframe.allocs <= Allocs(100)

    RNG = ComplexMixtures.init_random(options)
    t_RNG = @benchmark ComplexMixtures.init_random($options) samples = 1 evals = 1
    @test t_RNG.allocs <= Allocs(5)

    tmeta = ComplexMixtures.TrajectoryMetaData(traj, options)
    system = ComplexMixtures.ParticleSystem(traj, tmeta.unitcell, options, false, (1, 1))
    buff = ComplexMixtures.Buffer(traj, R)
    @. buff.solute_read = traj.x_solute
    @. buff.solvent_read = traj.x_solvent
    ComplexMixtures.update_unitcell!(system, ComplexMixtures.convert_unitcell(ComplexMixtures.getunitcell(traj)))
    t_mddf_frame =
        @benchmark ComplexMixtures.mddf_frame!($R, $system, $buff, $options, 1.0, $RNG) samples = 1 evals = 1
    @test t_mddf_frame.allocs < Allocs(100)

end
