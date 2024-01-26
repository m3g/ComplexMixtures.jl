"""

$(TYPEDEF)

Structure that contains the information about the solute and solvent molecules.

$(TYPEDFIELDS)

"""
@kwdef struct AtomSelection
    nmols::Int # Number of molecules
    natomspermol::Int # Number of atoms per molecule
    indices::Vector{Int} # Index of each atom in the full vector of coordinates
    # Groups of atoms: by default, the groups are the atoms themselves, but 
    # for very large molecules it is necessary to define groups of atoms 
    # in advance, to avoid memory issues. By default this is set to "nothing".
    custom_groups::Bool
    group_atom_indices::Vector{Vector{Int}}
    # Group (or atom) names.
    group_names::Vector{String} 
end

# AtomSelection show functions
function Base.show(io::IO, atsel::AtomSelection)
    (; nmols, natomspermol, custom_groups, group_atom_indices) = atsel 
    ngroups = custom_groups ? length(group_atom_indices) : natomspermol
    print(io, chomp(
    """
    AtomSelection 
        $(length(atsel.indices)) atoms belonging to $nmols molecule(s).
        Atoms per molecule: $natomspermol
        Number of groups: $ngroups
    """))
end

# Number of atoms of the selection
natoms(atsel::AtomSelection) = length(atsel.indices)

"""
    atom_group(atsel::AtomSelection, i::Int)
    atom_group(atsel::AtomSelection, groupname::String)

    atom_group(atsel::AtomSelection, i::Int)
    atom_group(atsel::AtomSelection, groupname::String)

Return the indices of the atoms that belong to a given group.

## Example

```jldoctest
julia> using ComplexMixtures

julia> atsel = AtomSelection([1,2,3], natomspermol=1, group_atom_indices=[[1,2],[3]], group_names=["G1", "G2"])
AtomSelection 
    3 atoms belonging to 3 molecule(s).
    Atoms per molecule: 1
    Number of groups: 2

julia> atom_group(atsel, 1)
2-element Vector{Int64}:
 1
 2

julia> atom_group(atsel, "G2")
1-element Vector{Int64}:
 3

julia> atom_group_name(atsel, 1)
"G1"
```

"""
atom_group(atsel::AtomSelection, i::Int) = atsel.group_atom_indices[i]
function atom_group(atsel::AtomSelection, group_name::String) 
    igroup = findfirst(==(group_name), atsel.group_names)
    if isnothing(igroup)
        error("Could not find group with name $group_name.")
    end
    return atsel.group_atom_indices[igroup]
end

"""
    atom_group_name(atsel::AtomSelection, i::Int)
    atom_group_names(atsel::AtomSelection)

Return the name of the group of atoms with index `i`. 
The `atom_group_names` function returns a vector with the names of all the groups.

## Example

```jldoctest
julia> using ComplexMixtures

julia> atsel = AtomSelection([1,2,3], natomspermol=1, group_atom_indices=[[1,2],[3]], group_names=["G1", "G2"])
AtomSelection 
    3 atoms belonging to 3 molecule(s).
    Atoms per molecule: 1
    Number of groups: 2

julia> atom_group_name(atsel, 1)
"G1"

julia> atom_group_names(atsel)
2-element Vector{String}:
 "G1"
 "G2"
```

"""
atom_group_name(atsel::AtomSelection, i::Int) = atsel.group_names[i]
@doc (@doc atom_group_name) atom_group_names(atsel) = atsel.group_names

"""

**AtomSelection constructors**

The `AtomSelection` structure carries the information of the molecules that are going to be used to compute the MDDF.
The structure can be initialized in different ways:

1. Initialize the structure providing a vector of PDBTools.Atom(s).

```
    AtomSelection(
        atoms::AbstractVector{<:PDBTools.Atom}; 
        nmols::Int = 0, 
        natomspermol::Int = 0,
        group_atom_indices::Union{Nothing,Vector{Vector{Int}}} = nothing,
        group_names::Vector{String} = String[]
    ) 
```

The indices of the atoms will be retrived from the
indices of the atoms as defined in the PDB file, thus the PDB file must correspond to the same
system as that of the simulation. 

Either the number of molecules (`nmols`) **or** the number of atoms per molecule (`natomspermol`) must be provided.

If `group_atom_indices` is `nothing` or `group_names` is empty, the names of the groups
will be retrieved from the atom names, and in the coordination numbers of each individual atom will be 
stored.

## Example

```jldoctest
julia> using ComplexMixtures, PDBTools

julia> pdbfile = ComplexMixtures.Testing.pdbfile;

julia> atoms = PDBTools.readPDB(pdbfile, "resname TMAO");

julia> atsel = AtomSelection(atoms, natomspermol=14)
AtomSelection 
    2534 atoms belonging to 181 molecule(s).
    Atoms per molecule: 14
    Number of groups: 14 

julia> atom_group_name(atsel, 1)
"N"

julia> atom_group_name(atsel, 5)
"O1"

julia> length(atom_group_names(atsel))
14
```

2. Lower level: initialize the structure providing the index of atoms and groups.

```
    AtomSelection(
        indices::Vector{Int};
        nmols::Int = 0,
        natomspermol::Int = 0,
        group_atom_indices::Union{Nothing,Vector{Vector{Int}}} = nothing,
        group_names::Vector{String} = String[]
    )
```

Construct an AtomSelection structure from the most low-level information: the index of atoms and groups.

Either the number of molecules (`nmols`) or the number of atoms per molecule (`natomspermol`) must be provided.

Groups of atoms can be defined by providing a vector of vectors of atom indices (`group_atom_indices`), and a vector of group names (`group_names`).
If `group_atom_indices` is set to `nothing`, the coordination numbers of each individual atoms wil be stored.

## Examples

```jldoctest
julia> using ComplexMixtures

julia> AtomSelection([1,2,3], nmols=1)
AtomSelection 
    3 atoms belonging to 1 molecule(s).
    Atoms per molecule: 3
    Number of groups: 3

julia> AtomSelection([1,2,3], natomspermol=1)
AtomSelection 
    3 atoms belonging to 3 molecule(s).
    Atoms per molecule: 1
    Number of groups: 1

julia> AtomSelection([1,2,3], natomspermol=1, group_atom_indices=[[1,2],[3]], group_names=["G1", "G2"])
AtomSelection 
    3 atoms belonging to 3 molecule(s).
    Atoms per molecule: 1
    Number of groups: 2 
```

"""
function AtomSelection(args...; kargs...) 
    throw(ArgumentError("No constructor for AtomSelection with these arguments. Please check the documentation."))
end

#
# Function to initialize the AtomSelection structure, from the 
# most low-level information: the index of atoms and groups
#
function AtomSelection(
    indices::Vector{Int};
    nmols::Int = 0,
    natomspermol::Int = 0,
    group_atom_indices::Vector{Vector{Int}} = Vector{Int}[],
    group_names::Vector{String} = String[]

)
    if nmols == 0 && natomspermol == 0
        throw(ArgumentError("Set nmols or natomspermol when defining a selection."))
    end
    natoms = length(indices)
    if natoms == 0
        throw(ArgumentError("Vector of atom indices is empty."))
    end
    if nmols != 0
        if natoms % nmols != 0
            throw(ArgumentError("Number of atoms in selection must be a multiple of nmols."))
        end
        natomspermol = div(natoms, nmols)
    else
        if natoms % natomspermol != 0
            throw(ArgumentError(" Number of atoms in selection must be a multiple of natomspermols."))
        end
        nmols = div(natoms, natomspermol)
    end

    # If the group_atom_indices is not empty, there are custom groups defined.
    custom_groups = !isempty(group_atom_indices)

    if !custom_groups && !isempty(group_names)
        if length(group_names) != natomspermol
            throw(ArgumentError(replace("""
            The length of the group_names vector does not correspond to the number of atoms per molecule,
            but no atom groups vector was provided. 
            """, '\n' => "")))
        end
    end

    if custom_groups
        if isempty(group_names)
            @warn begin
                """
                Vector of group atom indices was provided but vector of group names is empty. 
                The group contributions will be only retrieved by the group indices.
                """
            end _file=nothing _line=nothing
        else
            if length(group_atom_indices) != length(group_names)
                throw(ArgumentError("The vector of group atom indices has a different number of elements than the vector of group names."))
            end
        end
    end

    return AtomSelection(
        indices = indices,
        nmols = nmols,
        natomspermol = natomspermol,
        custom_groups = custom_groups,
        group_atom_indices = group_atom_indices,
        group_names = group_names,
    )
end

@testitem "AtomSelection - indices" begin
    using PDBTools: readPDB, index
    pdbfile = ComplexMixtures.Testing.pdbfile
    indices = index.(readPDB(pdbfile, "protein and residue 2"))
    s = AtomSelection(indices, nmols = 1, natomspermol = 11)
    @test s.indices == [12 + i for i = 1:11]
    @test s.group_names == String[]
    @test length(s.indices) == 11
    @test s.natomspermol == 11
    @test s.nmols == 1
    @test s.custom_groups == false
    @test s.group_names == String[]
    @test ComplexMixtures.natoms(s) == s.nmols * s.natomspermol
    s = AtomSelection(indices, group_names = fill("C", length(indices)), nmols = 1, natomspermol = 11)
    @test s.custom_groups == false
    @test s.group_names == fill("C", length(indices))
    @test ComplexMixtures.natoms(s) == s.nmols * s.natomspermol
end

@testitem "AtomSelection - argument errors" begin
    using ComplexMixtures
    @test_throws ArgumentError AtomSelection([1,2,3])
    @test_throws ArgumentError AtomSelection([1,2,3]; natomspermol=2)
    @test_throws ArgumentError AtomSelection([1,2,3]; natomspermol=1, nmols=2)
    @test_throws ArgumentError AtomSelection([1,2,3]; natomspermol=1, nmols=2)
    @test_throws ArgumentError AtomSelection([1,2,3]; natomspermol=1, group_names=["A", "B"])
    @test_throws ArgumentError AtomSelection([1,2,3]; natomspermol=1, group_names=["A", "B", "C"])
    
    @test_throws ArgumentError AtomSelection([1,2,3], ["A", "B", "C"])
    @test_throws MethodError AtomSelection([1,2,3]; abc = 1)
    @test_throws ArgumentError AtomSelection([1,2,3], ["A", "B", "C"]; nmols = 1)

    @test_logs (:warn,) AtomSelection([1,2,3], natomspermol=1, group_atom_indices=[[1,2],[3]])
end

#
# Initialize the structure providing a vector of PDBTools.Atom(s)
#
function AtomSelection(
    atoms::AbstractVector{<:PDBTools.Atom}; 
    nmols::Int = 0, 
    natomspermol::Int = 0,
    group_atom_indices::Vector{Vector{Int}} = Vector{Int}[],
    group_names::Vector{String} = String[]
)
    indices = PDBTools.index.(atoms)
    custom_groups = !isempty(group_atom_indices)
    if !custom_groups && isempty(group_names) 
        group_names = PDBTools.name.(atoms[1:natomspermol])
    end
    return AtomSelection(
        indices;
        nmols = nmols, 
        natomspermol = natomspermol,
        group_atom_indices = group_atom_indices,
        group_names = group_names, 
    )
end

@testitem "AtomSelection Vector{PDBTools.Atom}" begin
    import PDBTools
    pdbfile = ComplexMixtures.Testing.pdbfile
    atoms = PDBTools.readPDB(pdbfile, "protein and residue 2")
    s = AtomSelection(atoms, nmols = 1, natomspermol = 11)
    @test s.indices == [12 + i for i = 1:11]
    @test s.custom_groups == false
    @test atom_group_names(s) == ["N", "HN", "CA", "HA", "CB", "HB1", "HB2", "SG", "HG1", "C", "O"]
    @test length(s.indices) == 11
    @test s.natomspermol == 11
    @test s.nmols == 1
    @test ComplexMixtures.natoms(s) == s.nmols * s.natomspermol
end

"""
   **SoluteGroup** and **SolventGroup** data structures.

These structures are used to select groups of atoms to extract their contributions 
from the MDDF results. 

Most tipically, the groups are defined from a selection of atoms with the PDBTools package,
or by providing directly the indices of teh atoms in the structure. 

Alternativelly, if the groups were predefined, the groups can be selected by group index
or group name. 

The possible constructors are:

    SoluteGroup(atoms::Vector{PDBTools.Atom})
    SolventGroup(atoms::Vector{PDBTools.Atom})

    SoluteGroup(atom_indices::Vector{Int})
    SolventGroup(atom_indices::Vector{Int})

    SoluteGroup(atom_names::Vector{String})
    SolventGroup(atom_names::Vector{String})

    SoluteGroup(group_name::String) 
    SoluteGroup(group_index::Int) 

    SolventGroup(group_name::String)
    SolventGroup(group_index::Int)

    SolventGroup(residue::PDBTools.Residue)
    SolventGroup(residue::PDBTools.Residue)

""" SoluteGroup
@doc (@doc SoluteGroup) SolventGroup

struct SoluteGroup{
    I<:Union{Int,Nothing},
    S<:Union{String,Nothing},
    VI<:Union{Vector{Int},Nothing},
    VS<:Union{Vector{String},Nothing}
}
    group_index::I
    group_name::S
    atom_indices::VI
    atom_names::VS
end

struct SolventGroup{
    I<:Union{Int,Nothing},
    S<:Union{String,Nothing},
    VI<:Union{Vector{Int},Nothing},
    VS<:Union{Vector{String},Nothing}
}
    group_index::I
    group_name::S
    atom_indices::VI
    atom_names::VS
end

function SoluteGroup(args...; kargs...)
    throw(ArgumentError(chomp("""
        No constructor for SoluteGroup with these arguments. Please check the documentation.
    """)))
end
function SolventGroup(args...; kargs...)
    throw(ArgumentError(chomp("""
        No constructor for SoluteGroup with these arguments. Please check the documentation.
    """)))
end

SoluteGroup(atoms::Vector{PDBTools.Atom}) = SoluteGroup(nothing, nothing, PDBTools.index.(atoms), nothing)
SoluteGroup(atom_indices::Vector{Int}) = SoluteGroup(nothing, nothing, atom_indices, nothing)
SoluteGroup(atom_names::Vector{String}) = SoluteGroup(nothing, nothing, nothing, atom_names)
SoluteGroup(group_name::String) = SoluteGroup(nothing, group_name, nothing, nothing)
SoluteGroup(group_index::Int) = SoluteGroup(group_index, nothing, nothing, nothing)
SoluteGroup(residue::PDBTools.Residue) = SoluteGroup(nothing, nothing, PDBTools.index.(residue), nothing)

SolventGroup(atoms::Vector{PDBTools.Atom}) = SolventGroup(nothing, nothing, PDBTools.index.(atoms), nothing)
SolventGroup(atom_indices::Vector{Int}) = SolventGroup(nothing, nothing, atom_indices, nothing)
SolventGroup(atom_names::Vector{String}) = SolventGroup(nothing, nothing, nothing, atom_names)
SolventGroup(group_name::String) = SolventGroup(nothing, group_name, nothing, nothing)
SolventGroup(group_index::Int) = SolventGroup(group_index, nothing, nothing, nothing)
SolventGroup(residue::PDBTools.Residue) = SolventGroup(nothing, nothing, PDBTools.index.(residue), nothing)

@testitem "SoluteGroup and SolventGroup" begin
    using PDBTools: readPDB, select, name, eachresidue
    using ComplexMixtures.Testing: pdbfile
    pdb = readPDB(pdbfile)
    @test fieldnames(SoluteGroup) == fieldnames(SolventGroup)

    sg = SoluteGroup(select(pdb, "protein and residue 2"))
    @test sg.atom_indices == [12 + i for i = 1:11]
    @test count(!isnothing, getfield(sg, field) for field in fieldnames(SoluteGroup)) == 1
    @test SoluteGroup([1,2,3]).atom_indices == [1,2,3]
    @test count(!isnothing, getfield(sg, field) for field in fieldnames(SoluteGroup)) == 1
    @test SoluteGroup("N").group_name == "N"
    @test count(!isnothing, getfield(sg, field) for field in fieldnames(SoluteGroup)) == 1
    @test SoluteGroup(2).group_index == 2
    @test count(!isnothing, getfield(sg, field) for field in fieldnames(SoluteGroup)) == 1
    sg = SoluteGroup(name.(select(pdb, "protein and residue 2")))
    @test sg.atom_names == ["N", "HN", "CA", "HA", "CB", "HB1", "HB2", "SG", "HG1", "C", "O"]
    @test count(!isnothing, getfield(sg, field) for field in fieldnames(SoluteGroup)) == 1
    sg = SoluteGroup(collect(eachresidue(pdb))[2])
    @test sg.atom_indices == [12 + i for i = 1:11]
    @test count(!isnothing, getfield(sg, field) for field in fieldnames(SoluteGroup)) == 1

    sg = SolventGroup(select(pdb, "protein and residue 2"))
    @test sg.atom_indices == [12 + i for i = 1:11]
    @test count(!isnothing, getfield(sg, field) for field in fieldnames(SoluteGroup)) == 1
    @test SolventGroup([1,2,3]).atom_indices == [1,2,3]
    @test count(!isnothing, getfield(sg, field) for field in fieldnames(SoluteGroup)) == 1
    @test SolventGroup("N").group_name == "N"
    @test count(!isnothing, getfield(sg, field) for field in fieldnames(SoluteGroup)) == 1
    @test SolventGroup(2).group_index == 2
    @test count(!isnothing, getfield(sg, field) for field in fieldnames(SoluteGroup)) == 1
    sg = SolventGroup(name.(select(pdb, "protein and residue 2")))
    @test sg.atom_names == ["N", "HN", "CA", "HA", "CB", "HB1", "HB2", "SG", "HG1", "C", "O"]
    @test count(!isnothing, getfield(sg, field) for field in fieldnames(SoluteGroup)) == 1
    sg = SolventGroup(collect(eachresidue(pdb))[2])
    @test sg.atom_indices == [12 + i for i = 1:11]
    @test count(!isnothing, getfield(sg, field) for field in fieldnames(SoluteGroup)) == 1
end



