PrecompileTools.@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    dir = "$data_dir/NAMD"
    atoms = PDBTools.read_pdb("$dir/structure.pdb", "protein or resname TMAO")
    prot = PDBTools.select(atoms, "protein and resnum < 4")
    tmao = PDBTools.select(atoms, "resname TMAO and resnum <= 2")
    PrecompileTools.@compile_workload begin
        options = Options(lastframe=1, silent=true, n_random_samples=1)
        solute = AtomSelection(prot, nmols=1)
        solvent = AtomSelection(tmao, natomspermol=14)
        trajectory_file = "$dir/trajectory.dcd"
        R = mddf(trajectory_file, solute, solvent, options)
        rc = ResidueContributions(R, prot; silent=true)
        R = mddf(trajectory_file, solvent, options)
        tmpfile = tempname()
        save(tmpfile, R)
        load(tmpfile)
        rm(tmpfile)
    end
end
