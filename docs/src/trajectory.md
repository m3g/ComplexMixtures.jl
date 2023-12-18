# [Loading trajectories](@id trajectories)

To initialize a trajectory file for computation, use the command
```julia
trajectory = Trajectory("trajectory.xtc",solute,solvent)
```
where `solute` and `solvent` are defined with the `AtomSelection` function 
described [before](@ref selections). This function opens the stream for
reading frames, which are read once a time when the coordinates are
required for computing the MDDF.

The `Trajectory` function uses
[Chemfiles](http://chemfiles.org/Chemfiles.jl/latest/) in background,
and thus the most common trajectory formats are supported, as the ones
produced with NAMD, Gromacs, LAMMPS, Amber, etc.  

!!! tip
    The format of the trajectory file is automatically determined by
    Chemfiles from the extension of the file. However, it can be
    provided by the user with the `format` keyword, for example:
    ```julia
    trajectory = Trajectory("trajectory.xtc",solute,solvent,format="xtc")
    ```
