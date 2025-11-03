@testitem "Gromacs" begin
    using ComplexMixtures: mddf, Trajectory, Options, AtomSelection
    using ComplexMixtures: data_dir
    using PDBTools: read_pdb, select

    #
    # Tests with Gromacs-XSC trajectory
    #
    dir = "$data_dir/Gromacs"
    atoms = read_pdb("$dir/system.pdb")
    options = Options(stride=5, seed=321, StableRNG=true, nthreads=1, silent=true)

    # Example 1: protein-EMIM
    # save(R,"$dir/protein_EMI.json")
    R_save = load("$dir/protein_EMI.json")
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    emi = AtomSelection(select(atoms, "resname EMI"), natomspermol=20)
    traj = Trajectory("$dir/trajectory.xtc", protein, emi)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Example 2: EMIM-DCA
    # save(R,"$dir/EMI_DCA.json")
    R_save = load("$dir/EMI_DCA.json")
    emi = AtomSelection(select(atoms, "resname EMI"), natomspermol=20)
    dca = AtomSelection(select(atoms, "resname NC"), natomspermol=5)
    traj = Trajectory("$dir/trajectory.xtc", emi, dca)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

    # Example 3: EMIM-EMIM
    # save(R,"$dir/EMI_EMI.json")
    R_save = load("$dir/EMI_EMI.json")
    emi = AtomSelection(select(atoms, "resname EMI"), natomspermol=20)
    traj = Trajectory("$dir/trajectory.xtc", emi)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug=true)

end
