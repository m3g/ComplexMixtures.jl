
# [Set the solute and solvent selections](@id selections)

The solute and solvent are defined in ComplexMixtures as lists (vectors) of the
indexes of the atoms of the system. The solute and solvent information
is stored in the `Selection` structure. For example, if the solute is a
molecule formed by the first 5 atoms of the system, it would be defined
as:     
```julia
indexes = [ 1, 2, 3, 4, 5 ]
solute = Selection(indexes,nmols=1)
```

!!! note
    We need to inform the `Selection` function about the number of atoms of
    each molecule (using `natomspermol=3`, for example), or the number 
    of molecules (using `nmols=1000`, for example),
    such that the atoms belonging to each molecule can be determined
    without ambiguity. 

The atom names can be also provided such that some of the output files
contain more information on the [atomic contributions](@ref contrib). In this
case the syntax is:
```julia
indexes = [ 1, 2, 3, 4, 5 ]
names = [ "H1", "H2", "H3", "H4", "C" ]
solute = Selection(indexes,names,nmols=1)
```

!!! warning
    The indexing in ComplexMixtures is 1-based. That means that the first atom of
    your structure file is in position 1 of the coordinates. Please be
    careful if using any selection tool to be sure that your selection
    is correct.


## Using PDBTools

[PDBTools](https://github.com/m3g/PDBTools) is a package we developed to read and 
write PDB files,
which provides a simple selection tool. It is installed as a dependency 
of ComplexMixtures.  Given a PDB file of the simulated system, the solute can
be defined as, for example,
```julia
using PDBTools
atoms = PDBTools.readPDB("system.pdb")
protein = PDBTools.select(atoms,"protein")
solute = Selection(protein,nmols=1)
```
If the solvent is, for instance, water, the indexes of the water
molecules can be obtained with:
```julia
water = PDBTools.select(atoms,"water")
solvent = Selection(water,natomspermol=3)
```
or, alternatively, a more compact syntax can be used, for example:
```julia
water = PDBTools.select("system.pdb","resname TIP3P")
solvent = Selection(water,natomspermol=3)
```

or even providing just the names of the input file and selection, which
will run PDBTools in background:
```julia
solvent = Selection("sytem.pdb","water",water,natomspermol=3)
```
## Using VMD

[VMD](https://www.ks.uiuc.edu/Research/vmd/) is a very popular and
powerful package for visualization of simulations. It contains a very
versatile library to read toppologies and trajectory files, and a
powerful selection syntax. We provide here a wrapper to VMD which allows
using its capabilities.  

For example, the solute can be defined with: 
```julia
indexes, names = VMDselect("./system.gro","protein",vmd="/usr/bin/vmd")
solute = Selection(indexes,names,nmols=1)
```
The main advantage here is that all the file types that VMD supports are
supported. But VMD needs to be installed and is run in background, and
it takes a few seconds.     

!!! warning
    VMD uses 0-based indexing and `VMDselect` adjusts that. However, if
    a selection is performed by index, as with `index 1`, VMD will
    select the second atom, and the output will be `[2]`. Selections by
    type, name, segment, residue name, etc, won't be a problem.

