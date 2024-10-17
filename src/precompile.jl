PrecompileTools.@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    dir = "$(Testing.data_dir)/NAMD"
    atoms = PDBTools.readPDB("$dir/structure.pdb", "protein or resname TMAO")
    prot = PDBTools.select(atoms, "protein and resnum < 4")
    tmao = PDBTools.select(atoms, "resname TMAO and resnum <= 2")
    PrecompileTools.@compile_workload begin
        options = Options(lastframe = 1, silent = true, n_random_samples=1)
        solute = AtomSelection(prot, nmols = 1)
        solvent = AtomSelection(tmao, natomspermol = 14)
        traj = Trajectory("$dir/trajectory.dcd", solute, solvent)
        R = mddf(traj, options)
        rc = ResidueContributions(R, prot; silent=true)
        traj = Trajectory("$dir/trajectory.dcd", solvent)
        R = mddf(traj, options)
    end
end
