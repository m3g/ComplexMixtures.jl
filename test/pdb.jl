
@testitem "PDB" begin
    using ComplexMixtures: AtomSelection, Options, mddf, Trajectory, load
    using PDBTools: readPDB, select
    using ComplexMixtures.Testing: data_dir
    #
    # Tests with trajectory given in a PDB file
    #
    dir = "$data_dir/PDB"
    atoms = readPDB("$dir/trajectory.pdb", "model 1")
    options = Options(stride=1, seed=321, StableRNG=true, nthreads=1, silent=true)

    # Example 1: protein-tmao
    # save(R,"$dir/protein_tmao.json")
    R_save = load("$dir/protein_tmao.json")
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    traj = Trajectory("$dir/trajectory.pdb", protein, tmao, format="PDBTraj")
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Example 2: tmao-tmao
    # save(R,"$dir/tmao_tmao.json")
    R_save = load("$dir/tmao_tmao.json")
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    traj = Trajectory("$dir/trajectory.pdb", tmao, format="PDBTraj")
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)
end
