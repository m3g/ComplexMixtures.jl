
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

julia> atoms = read_pdb(ComplexMixtures.Testing.pdbfile);

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

```julia-repl
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
powerful selection syntax. The PDBTools.jl (v1.0 or greater) package provides a simple
wrapper to VMD that allows using the same syntax at it supports.

For example, the solute can be defined with: 
```julia
using ComplexMixtures, PDBTools

atoms = read_pdb("system.pdb")

indices, names = select_with_vmd("./system.pdb", "protein", vmd="/usr/bin/vmd")

protein = atoms[indices]

solute = AtomSelection(protein, nmols=1)
```
The main advantage here is that all the file types that VMD supports are
supported. But VMD needs to be installed and is run in background, and
it takes a few seconds to be executed. 

The `VMDSelect` function also accepts an optional keyword parameter `srcload`,
which can be used to load custom scripts within `vmd` before setting
the selection. This allows the definition of `tcl` scripts with custom selection
macros, for instance. The usage would be: 
```julia
using PDBTools

indices, names = select_with_vmd(
    "file.pdb", 
    "resname MYRES"; 
    srcload = [ "mymacros1.tcl", "mymacros2.tcl" ]
)
```
Which corresponds to `source`ing each of the macro files in VMD before defining the 
selection with the custom `MYRES` name.

!!! warning
    VMD uses 0-based indexing and `VMDselect` adjusts that. However, if
    a selection is performed by index, as with `index 1`, VMD will
    select the second atom, and the output will be `[2]`. AtomSelections by
    type, name, segment, residue name, etc, won't be a problem.

## [Predefinition of atom groups](@id predefinition-of-groups)

Importantly, this should be only a concern for the solvation analysis of systems in which individual
molecules are **very large**. This feature was introduced in version `2.0` of the package to support
the study of small molecule distribution in virus structures, of millions of atoms.  

By default, the contribution of each *type of* atom to the coordination number counts is stored, to allow
the decomposition of the final MDDFs into any group contribution. However, when a structure, like a virus,
has millions of atoms, storing the contribution of each atom becomes prohibitive in terms of memory.
Thus, one may need to predefine the groups in which the contributions will be analyzed.

Here, we illustrate this feature by presselecting the acidic and basic residues of a protein:

```julia
julia> using ComplexMixtures, PDBTools

julia> atoms = read_pdb(ComplexMixtures.Testing.pdbfile);

julia> protein = select(atoms, "protein");

julia> acidic_residues = select(atoms, "protein and acidic");

julia> basic_residues = select(atoms, "protein and basic");

julia> solute = AtomSelection(
        protein, 
        nmols=1,
        group_atom_indices = [ index.(acidic_residues), index.(basic_residues) ],
        group_names = [ "acidic residues", "basic residues" ]
    )
AtomSelection 
    1463 atoms belonging to 1 molecule(s).
    Atoms per molecule: 1463
    Number of groups: 2 
```

In this example, then, the `solute` `AtomSelection` has two groups. The indices of the atoms
of the groups are stored in the `group_atom_indices` vector and the group names in the `group_names`
vector. The `atom_group` auxiliary function is the most practical way to retrive the indices of the
atoms of the group.

```julia-repl
julia> atom_group(solute, "acidic residues")
162-element Vector{Int64}:
   24
   25
   26
    ⋮
 1436
 1437
```

With these group selections predefined, the contributions of these groups to the MDDF or coordination numbers
can be retrived directly from the result data structure with, for example:

```julia-repl
julia> result = mddf(trajectory_file, solute, solvent, Options(bulk_range=(8.0, 12.0)));

julia> acidic_residue_contributions = contributions(result, SoluteGroup("acidic residues"))
```

## Reference functions

```@autodocs
Modules = [ComplexMixtures]
Pages = ["AtomSelection.jl"]
```








