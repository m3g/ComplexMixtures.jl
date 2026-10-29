
@testitem "NAMD" begin
    using ComplexMixtures
    using LazyArtifacts
    using ComplexMixtures: test_file
    using PDBTools: read_pdb, select, Select

    #
    # Tests with NAMD-DCD trajectory
    #
    atoms = read_pdb(test_file("protein_tmao.pdb"))
    options = Options(stride=5, seed=321, StableRNG=true, nthreads=1, silent=true, cutoff=10.0)

    # Example 1: protein-tmao
    # save(R,test_file("namd_protein_tmao.json"))
    R_save = load(MDDFResult, test_file("namd_protein_tmao.json"))
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    traj = ComplexMixtures.Trajectory(test_file("protein_tmao.dcd"), protein, tmao)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # test substring input
    traj = ComplexMixtures.Trajectory(@view((test_file("protein_tmao.dcd")*"---")[1:end-3]), protein, tmao)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Test save and load
    temp_output = tempname()
    save(R, temp_output)
    R_load = load(MDDFResult, temp_output)
    @test R_load ≈ R_save
    rm(temp_output)

    # Example 2: water-tmao
    # save(R,test_file("namd_water_tmao.json"))
    R_save = load(MDDFResult, test_file("namd_water_tmao.json"))
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    water = AtomSelection(select(atoms, "water"), natomspermol=3)
    traj = ComplexMixtures.Trajectory(test_file("protein_tmao.dcd"), tmao, water)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Example 3: tmao-tmao
    # save(R,test_file("namd_tmao_tmao.json"))
    R_save = load(MDDFResult, test_file("namd_tmao_tmao.json"))
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    traj = ComplexMixtures.Trajectory(test_file("protein_tmao.dcd"), tmao)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Example 3: water-water
    # save(R,test_file("namd_water_water.json"))
    R_save = load(MDDFResult, test_file("namd_water_water.json"))
    water = AtomSelection(select(atoms, "water"), natomspermol=3)
    traj = ComplexMixtures.Trajectory(test_file("protein_tmao.dcd"), water)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Test show_progress path
    traj = ComplexMixtures.Trajectory(test_file("protein_tmao.dcd"), protein, tmao; show_progress=true)
    @test traj.nframes == 20
    traj = ComplexMixtures.Trajectory(test_file("protein_tmao.dcd"), protein, tmao; show_progress=false)
    @test traj.nframes == 20
end
