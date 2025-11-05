
@testitem "NAMD with ChemFiles" begin
    using ComplexMixtures
    using ComplexMixtures: data_dir
    using PDBTools: read_pdb, select

    #
    # Tests with NAMD-DCD trajectory
    #

    dir = "$data_dir/NAMD"
    atoms = read_pdb("$dir/structure.pdb")
    options = Options(stride=5, seed=321, StableRNG=true, nthreads=1, silent=true, bulk_range=(8.0, 10.0))

    # Example 1: protein-tmao
    # save(R,"$dir/protein_tmao.json")
    R_save = load("$dir/protein_tmao.json")
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    traj = Trajectory("$dir/trajectory.dcd", protein, tmao, chemfiles=true)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Test save and load
    temp_output = tempname()
    save(R, temp_output)
    R_load = load(temp_output)
    @test R_load ≈ R_save

    # Save and load a coordination-number run
    C = coordination_number(traj)
    save(C, temp_output)
    C_read = load(temp_output)
    @test C_read ≈ C
    rm(temp_output)

    # Example 2: water-tmao
    # save(R,"$dir/water_tmao.json")
    R_save = load("$dir/water_tmao.json")
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    water = AtomSelection(select(atoms, "water"), natomspermol=3)
    traj = Trajectory("$dir/trajectory.dcd", tmao, water, chemfiles=true)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Example 3: tmao-tmao
    # save(R,"$dir/tmao_tmao.json")
    R_save = load("$dir/tmao_tmao.json")
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    traj = Trajectory("$dir/trajectory.dcd", tmao, tmao, chemfiles=true)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Example 3: water-water
    # save(R,"$dir/water_water.json")
    R_save = load("$dir/water_water.json")
    water = AtomSelection(select(atoms, "water"), natomspermol=3)
    traj = Trajectory("$dir/trajectory.dcd", water, water, chemfiles=true)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

end
