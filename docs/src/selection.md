
# [Set the solute and solvent selections](@id selections)

The solute and solvent are defined in MDDF as lists (vectors) of the
indexes of the atoms of the system. The solute and solvent information
is stored in the `Selection` structure. For example, if the solute is a
molecule formed by the first 5 atoms of the system, it would be defined
as:     
```julia
indexes = [ 1, 2, 3, 4, 5 ]
solute = MDDF.Selection(indexes,nmols=1)
```
We need to inform, as well, how many molecules compose the solute, in
this case 1. Of course, writing explicit the lists of atoms can be
cumbersome, and selection tools are provided.      

!!! warning
    The indexing in MDDF is 1-based. That means that the first atom of
    your structure file is in position 1 of the coordinates. Please be
    careful if using any selection tool to be sure that your selection
    is correct.


## Using PDBTools

[PDBTools](https://github.com/m3g/PDBTools) is a simple package we developed to read and write PDB files,
which provides a simple selection tool. Install it according to the
instructions. Given a PDB file of the simulated system, the solute can
be defined as, for example,
```julia
indexes = PDBTools.select("system.pdb","protein")
solute = MDDF.Selection(indexes,nmols=1)
```
If the solvent is, for instance, water, the indexes of the water
molecules can be obtained with:
```julia
indexes = PDBTools.select("system.pdb","water")
solvent = MDDF.Selection(indexes,natomspermol=3)
```
or, alternatively,
```julia
indexes = PDBTools.select("system.pdb","resname TIP3P")
solvent = MDDF.Selection(indexes,natomspermol=3)
```
We need to inform the `Selection` function about the number of atoms of
each molecule, or the number of molecules (one information is obtained
from the other), such that the identity of the atoms of each molecule
can be obtained.   

## Using VMD

[VMD](https://www.ks.uiuc.edu/Research/vmd/) is a very popular and
powerful package for visualization of simulations. It contains a very
versatile library to read topologies and trajectory files, and a
powerful selection syntax. We provide here a wrapper to VMD which allows
using its capabilities.  

For example, the solute can be defined with: 
```julia
indexes = MDDF.VMDSelect("./system.gro","protein",vmd="/usr/bin/vmd")
solute = MDDF.Selection(indexes,nmols=1)
```
The main advantage here is that all the file types that VMD supports are
supported. But VMD needs to be installed and is run in background, and
it takes a few seconds.     

!!! warning
    VMD uses 0-based indexing and `VMDSelect` adjusts that. However, if
    a selection is performed by index, as with `index 1`, VMD will
    select the second atom, and the output will be `[2]`. Selections by
    type, name, segment, residue name, etc, won't be a problem.

