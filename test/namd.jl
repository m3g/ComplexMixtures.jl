
@testitem "NAMD" begin
    using ComplexMixtures: mddf, Trajectory, Options, AtomSelection, load,
        SoluteGroup, SolventGroup, contributions
    using PDBTools: readPDB, select, Select
    using ComplexMixtures.Testing: data_dir

    #
    # Tests with NAMD-DCD trajectory
    #
    dir = "$data_dir/NAMD"
    atoms = readPDB("$dir/structure.pdb")
    options = Options(stride=5, seed=321, StableRNG=true, nthreads=1, silent=true, bulk_range=(8.0, 10.0))

    # Example 1: protein-tmao
    # save(R,"$dir/protein_tmao.json")
    R_save = load("$dir/protein_tmao.json")
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    traj = Trajectory("$dir/trajectory.dcd", protein, tmao)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Test save and load
    temp_output = tempname()
    save(R, temp_output)
    R_load = load(temp_output)
    @test R_load ≈ R_save

    # Example 2: water-tmao
    # save(R,"$dir/water_tmao.json")
    R_save = load("$dir/water_tmao.json")
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    water = AtomSelection(select(atoms, "water"), natomspermol=3)
    traj = Trajectory("$dir/trajectory.dcd", tmao, water)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Example 3: tmao-tmao
    # save(R,"$dir/tmao_tmao.json")
    R_save = load("$dir/tmao_tmao.json")
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    traj = Trajectory("$dir/trajectory.dcd", tmao)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Example 3: water-water
    # save(R,"$dir/water_water.json")
    R_save = load("$dir/water_water.json")
    water = AtomSelection(select(atoms, "water"), natomspermol=3)
    traj = Trajectory("$dir/trajectory.dcd", water)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Test show_progress path
    traj = Trajectory("$dir/trajectory.dcd", protein, tmao; show_progress=true)
    @test traj.nframes == 20
    traj = Trajectory("$dir/trajectory.dcd", protein, tmao; show_progress=false)
    @test traj.nframes == 20
end
