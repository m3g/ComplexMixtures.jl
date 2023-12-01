@testitem "Gromacs" begin
    using ComplexMixtures, PDBTools
    using ComplexMixtures.Testing
    const CM = ComplexMixtures

    #
    # Tests with Gromacs-XSC trajectory
    #

    dir = "$(Testing.data_dir)/Gromacs"
    atoms = readPDB("$dir/system.pdb")
    options = Options(stride = 5, seed = 321, StableRNG = true, nthreads = 1, silent = true)

    # Example 1: protein-EMIM

    # save(R,"$dir/protein_EMI.json")
    R_save = load("$dir/protein_EMI.json"; legacy_warning = false)
    protein = Selection(select(atoms, "protein"), nmols = 1)
    emi = Selection(select(atoms, "resname EMI"), natomspermol = 20)
    traj = Trajectory("$dir/trajectory.xtc", protein, emi)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug = true)

    # Example 1: EMIM-DCA

    # save(R,"$dir/EMI_DCA.json")
    R_save = load("$dir/EMI_DCA.json"; legacy_warning = false)
    emi = Selection(select(atoms, "resname EMI"), natomspermol = 20)
    dca = Selection(select(atoms, "resname NC"), natomspermol = 5)
    traj = Trajectory("$dir/trajectory.xtc", emi, dca)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug = true)

    # Example 1: EMIM-EMIM

    # save(R,"$dir/EMI_EMI.json")
    R_save = load("$dir/EMI_EMI.json"; legacy_warning = false)
    emi = Selection(select(atoms, "resname EMI"), natomspermol = 20)
    traj = Trajectory("$dir/trajectory.xtc", emi)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug = true)

end
