
@testitem "PDB" begin
    using ComplexMixtures: AtomSelection, Options, mddf, Trajectory, load
    using LazyArtifacts
    using ComplexMixtures: test_file
    using PDBTools: read_pdb, select
    #
    # Tests with trajectory given in a PDB file
    #
    atoms = read_pdb(test_file("pdb_trajectory.pdb"), "model 1")
    options = Options(stride=1, seed=321, StableRNG=true, nthreads=1, silent=true)

    # Example 1: protein-tmao
    # save(R,test_file("pdb_protein_tmao.json"))
    R_save = load(MDDFResult, test_file("pdb_protein_tmao.json"))
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    traj = ComplexMixtures.Trajectory(test_file("pdb_trajectory.pdb"), protein, tmao, format="PDBTraj")
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Example 2: tmao-tmao
    # save(R,test_file("pdb_tmao_tmao.json"))
    R_save = load(MDDFResult, test_file("pdb_tmao_tmao.json"))
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    traj = ComplexMixtures.Trajectory(test_file("pdb_trajectory.pdb"), tmao, format="PDBTraj")
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)
end
