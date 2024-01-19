
# [Solute and solvent selections](@id selections)

The solute and solvent are defined by selecting subsets of atoms from the 
system. These subsets are defined by the `AtomSelection` data structures. 

To construct a `AtomSelection` data structure, one needs to provide, at least,
the (1-based) indices of the atoms that belong to the selection, and either
the number of atoms of each molecule or the number of molecules in the selection.

## Using the PDBTools package

The [PDBTools](https://m3g.github.io/PDBTools) package helps the construction of 
the solute and solvent data structures,
by providing a convenient selection syntax. Additionally, it sets up the names
of the atoms of the system in the data structure, which can be used to retrieve
atom and and group contributions to MDDFs and coordination numbers. 

For example, here we define a protein of a system as the solute:

```jldoctest
julia> using ComplexMixtures, PDBTools

julia> atoms = readPDB(ComplexMixtures.Testing.pdbfile);

julia> protein = select(atoms, "protein");

julia> solute = AtomSelection(protein, nmols=1)
AtomSelection 
    1463 atoms belonging to 1 molecule(s).
    Atoms per molecule: 1463
    Number of groups: 1463 
```

We need to inform the `AtomSelection` function about the number of atoms of
each molecule (using `natomspermol=3`, for example), or the number 
of molecules (using `nmols=1000`, for example),
such that the atoms belonging to each molecule can be determined
without ambiguity. 

Now, we define the solvent of the system as the water molecules:

```jldoctest
julia> water = select(atoms, "water"); 

julia> solvent = AtomSelection(water, natomspermol=3)
AtomSelection 
    58014 atoms belonging to 19338 molecule(s).
    Atoms per molecule: 3
    Number of groups: 3
```

## Using VMD

[VMD](https://www.ks.uiuc.edu/Research/vmd/) is a very popular and
powerful package for visualization of simulations. It contains a very
versatile library to read topologies and trajectory files, and a
powerful selection syntax. We provide here a wrapper to VMD which allows
using its capabilities.  

For example, the solute can be defined with: 
```julia
indexes, names = VMDselect("./system.gro","protein",vmd="/usr/bin/vmd")
solute = AtomSelection(indexes,names,nmols=1)
```
The main advantage here is that all the file types that VMD supports are
supported. But VMD needs to be installed and is run in background, and
it takes a few seconds to be excuted.

The `VMDSelect` function also accepts an optional keyword parameter `srcload`,
which can be used to load custom scripts within `vmd` before setting
the selection. This allows the definition of `tcl` scripts with custom selection
macros, for instance. The usage would be: 
```julia
sel = VMSelect("file.pdb", "resname MYRES"; srcload = [ "mymacros1.tcl", "mymacros2.tcl" ])
```
Which corresponds to `source`ing each of the macro files in VMD before defining the 
selection with the custom `MYRES` name.

!!! compat
    Custom script source loading in VMDSelect was introduced in ComplexMixtures version 1.3.0.

!!! warning
    VMD uses 0-based indexing and `VMDselect` adjusts that. However, if
    a selection is performed by index, as with `index 1`, VMD will
    select the second atom, and the output will be `[2]`. AtomSelections by
    type, name, segment, residue name, etc, won't be a problem.

## Predefinition of atom groups




