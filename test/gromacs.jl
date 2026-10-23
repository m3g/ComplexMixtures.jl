@testitem "Gromacs" begin
    using ComplexMixtures
    using LazyArtifacts
    using ComplexMixtures: test_file
    using PDBTools: read_pdb, select

    #
    # Tests with Gromacs-XSC trajectory
    #
    atoms = read_pdb(test_file("gromacs_system.pdb"))
    options = Options(stride=5, seed=321, StableRNG=true, nthreads=1, silent=true)

    # Example 1: protein-EMIM
    # save(R,test_file("gromacs_protein_EMI.json"))
    R_save = load(test_file("gromacs_protein_EMI.json"))
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    emi = AtomSelection(select(atoms, "resname EMI"), natomspermol=20)
    traj = Trajectory(test_file("gromacs_trajectory.xtc"), protein, emi)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # test substring input
    traj = Trajectory(@view((test_file("gromacs_trajectory.xtc")*"---")[1:end-3]), protein, emi)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Example 2: EMIM-DCA
    # save(R,test_file("gromacs_EMI_DCA.json"))
    R_save = load(test_file("gromacs_EMI_DCA.json"))
    emi = AtomSelection(select(atoms, "resname EMI"), natomspermol=20)
    dca = AtomSelection(select(atoms, "resname NC"), natomspermol=5)
    traj = Trajectory(test_file("gromacs_trajectory.xtc"), emi, dca)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Example 3: EMIM-EMIM
    # save(R,test_file("gromacs_EMI_EMI.json"))
    R_save = load(test_file("gromacs_EMI_EMI.json"))
    emi = AtomSelection(select(atoms, "resname EMI"), natomspermol=20)
    traj = Trajectory(test_file("gromacs_trajectory.xtc"), emi)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

end
