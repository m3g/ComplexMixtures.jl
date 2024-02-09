
@testitem "NAMD" begin
    using ComplexMixtures: mddf, Trajectory, Options, AtomSelection, load
    using PDBTools: readPDB, select
    using ComplexMixtures.Testing: data_dir

    #
    # Tests with NAMD-DCD trajectory
    #

    dir = "$data_dir/NAMD"
    atoms = readPDB("$dir/structure.pdb")
    options = Options(stride = 5, seed = 321, StableRNG = true, nthreads = 1, silent = true)

    # Example 1: protein-tmao
    # save(R,"$dir/protein_tmao.json")
    R_save = load("$dir/protein_tmao.json"; legacy_warning = false)
    protein = AtomSelection(select(atoms, "protein"), nmols = 1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol = 14)
    traj = Trajectory("$dir/trajectory.dcd", protein, tmao)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug = true)

    # Test save and load
    temp_output = tempname()
    save(R, temp_output)
    R_load = load(temp_output)
    @test R_load ≈ R_save

    # Example 2: water-tmao
    # save(R,"$dir/water_tmao.json")
    R_save = load("$dir/water_tmao.json"; legacy_warning = false)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol = 14)
    water = AtomSelection(select(atoms, "water"), natomspermol = 3)
    traj = Trajectory("$dir/trajectory.dcd", tmao, water)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug = true)

    # Example 3: tmao-tmao
    # save(R,"$dir/tmao_tmao.json")
    R_save = load("$dir/tmao_tmao.json"; legacy_warning = false)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol = 14)
    traj = Trajectory("$dir/trajectory.dcd", tmao)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug = true)

    # Example 3: water-water
    # save(R,"$dir/water_water.json")
    R_save = load("$dir/water_water.json"; legacy_warning = false)
    water = AtomSelection(select(atoms, "water"), natomspermol = 3)
    traj = Trajectory("$dir/trajectory.dcd", water)
    R = mddf(traj, options)
    @test isapprox(R, R_save, debug = true)

    #
    # Test computation of custom group counters
    #
    solute = AtomSelection(
        select(atoms , "water and residue 301"), 
        nmols = 1,
        group_atom_indices = [
            findall(Select("water and residue 301 and name H1"), atoms),
            findall(Select("water and residue 301 and name H2"), atoms),
            findall(Select("water and residue 301 and name OH2"), atoms)
        ],
        group_names = ["H1", "H2", "OH2"] 
    )
    solvent = AtomSelection(
        select(atoms, "water and not residue 301"),
        natomspermol = 3, 
        group_atom_indices = [
            findall(Select("water and name H1 and not residue 301"), atoms),
            findall(Select("water and name H2 and not residue 301"), atoms),
            findall(Select("water and name OH2 and not residue 301"), atoms)
        ],
        group_names = ["H1", "H2", "OH2"] 
    )
    traj = Trajectory("$dir/trajectory.dcd", solute, solvent)
    R = mddf(traj)

    @test !all(==(0), contributions(R, SoluteGroup("H1")))
    @test !all(==(0), contributions(R, SoluteGroup("H2")))
    @test !all(==(0), contributions(R, SoluteGroup("OH2")))

    @test !all(==(0), contributions(R, SolventGroup("H1")))
    @test !all(==(0), contributions(R, SolventGroup("H2")))
    @test !all(==(0), contributions(R, SolventGroup("OH2")))

    @test contributions(R, SoluteGroup("H1")) + 
          contributions(R, SoluteGroup("H2")) + 
          contributions(R, SoluteGroup("OH2")) ≈ 
          contributions(R, SolventGroup("H1")) + 
          contributions(R, SolventGroup("H2")) + 
          contributions(R, SolventGroup("OH2"))

    @test contributions(R, SoluteGroup("H1")) + 
          contributions(R, SoluteGroup("H2")) + 
          contributions(R, SoluteGroup("OH2")) ≈ 
          R.mddf

end
