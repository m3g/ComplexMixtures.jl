
# [Loading trajectories](@id trajectories)

To initialize a trajectory file for computation, use the command
```julia
trajectory = ComplexMixtures.Trajectory("trajectory.xtc",solute,solvent)
```
where `solute` and `solvent` are defined with the `Selection` function 
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
    trajectory = ComplexMixtures.Trajectory("trajectory.xtc",solute,solvent,format="xtc")
    ```

!!! note 
    The trajectory stream is closed at the end of the MDDF computation.
    Therefore if you want to reuse the same trajectory for another MDDF 
    computation in the same script, you need to reload it. For example:
    ```julia
    solute = ComplexMixtures.Selection("system.pdb","protein",nmols=1)
    # Compute the protein-water MDDF
    solvent = ComplexMixtures.Selection("system.pdb","water",natomspermol=3)
    trajectory = ComplexMixtures.Trajectory("trajectory.xtc",solute,solvent)
    R_water = ComplexMixtures.mddf(trajectory)
    # Compute the protein-urea MDDF
    solvent = ComplexMixtures.Selection("system.pdb","resname URE",natomspermol=8)
    trajectory = ComplexMixtures.Trajectory("trajectory.xtc",solute,solvent)
    R_urea = ComplexMixtures.mddf(trajectory)
    ```



