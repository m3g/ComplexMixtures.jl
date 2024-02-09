PrecompileTools.@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    dir = "$(Testing.data_dir)/NAMD"
    atoms = PDBTools.readPDB("$dir/structure.pdb", "resname TMAO")
    tmao1 = AtomSelection(PDBTools.select(atoms, "resname TMAO and resnum 1"), natomspermol = 14)
    PrecompileTools.@compile_workload begin
        options = Options(lastframe = -1, silent = true)
        tmao2 = AtomSelection(
            PDBTools.select(atoms, "resname TMAO and resnum 2 or resname TMAO and resnum 3"),
            natomspermol = 14,
        )
        traj = Trajectory("$dir/trajectory.dcd", tmao1, tmao2)
        R = mddf(traj, options)
        traj = Trajectory("$dir/trajectory.dcd", tmao2)
        R = mddf(traj, options)
    end
end
