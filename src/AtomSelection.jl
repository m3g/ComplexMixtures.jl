"""

$(TYPEDEF)

Structure that contains the information about the solute and solvent molecules.

$(TYPEDFIELDS)

"""
@kwdef struct AtomSelection
    natoms::Int # Total number of atoms
    nmols::Int # Number of molecules
    natomspermol::Int # Number of atoms per molecule
    index::Vector{Int} # Index of each atom in the full vector of coordinates
    imol::Vector{Int} # index of the molecule to which each atom belongs
    # Groups of atoms: by default, the groups are the atoms themselves, but 
    # for very large molecules it is necessary to define groups of atoms 
    # in advance, to avoid memory issues. By default this vector is empty.
    group::Vector{Vector{Int}}
    # Group (or atom) names.
    name::Vector{String} 
end

# Initialize providing the file name, and calling by default PDBTools.select
function AtomSelection(
    pdbfile::String, 
    selection::String; 
    nmols::Int = 0, 
    natomspermol::Int = 0,
    groups::AbstractVector{AbstractVector{Int}} = Vector{Int}[],
    names::AbstractVector{String} = String[]
)
    sel = PDBTools.readPDB(pdbfile, selection)
    return AtomSelection(
        sel, 
        nmols = nmols, 
        natomspermol = natomspermol,
        names = names,
        groups = groups,
    )
end

@testitem "AtomSelection PDBTools" begin
    pdbfile = ComplexMixtures.Testing.pdbfile
    s = AtomSelection(pdbfile, "protein and residue 2", nmols = 1, natomspermol = 11)
    @test s.imol == ones(Int, 11)
    @test s.index == [12 + i for i = 1:11]
    @test s.names == ["N", "HN", "CA", "HA", "CB", "HB1", "HB2", "SG", "HG1", "C", "O"]
    @test s.natoms == 11
    @test s.natomspermol == 11
    @test s.nmols == 1
end

# If the input is a vector of PDBTools.Atom types, load the index and types
function AtomSelection(
    atoms::AbstractVector{<:PDBTools.Atom}; 
    nmols::Int = 0, 
    natomspermol::Int = 0,
    groups::Vector{Vector{Int}} = Vector{Int}[],
    names::Vector{String} = String[]
)
    index = PDBTools.index.(atoms)
    if isempty(groups) && !isempty(names)
        throw(ArgumentError("Vector of group names was provided but vector of groups is empty."))
    end 
    if !isempty(groups) && isempty(names)
        @warn begin
            """
            Vector of groups was provided but vector of group names is empty. 
            The group contributions will be only retrieved by the group index.
            """
        end _file=nothing _line=nothing
    end 
    if isempty(groups)
        name = PDBTools.name.(atoms)
    end
    return AtomSelection(
        natoms = length(atoms),
        index = index,
        name = name, 
        nmols = nmols, 
        natomspermol = natomspermol,
        groups = groups
    )
end

@testitem "AtomSelection Vector{PDBTools.Atom}" begin
    import PDBTools
    pdbfile = ComplexMixtures.Testing.pdbfile
    atoms = PDBTools.readPDB(pdbfile, "protein and residue 2")
    s = AtomSelection(atoms, nmols = 1, natomspermol = 11)
    @test s.imol == ones(Int, 11)
    @test s.index == [12 + i for i = 1:11]
    @test s.names == ["N", "HN", "CA", "HA", "CB", "HB1", "HB2", "SG", "HG1", "C", "O"]
    @test s.natoms == 11
    @test s.natomspermol == 11
    @test s.nmols == 1
end

# If no names are provided, just repeat the indexes
function AtomSelection(indexes::Vector{Int}; nmols::Int = 0, natomspermol::Int = 0)
    names = ["$i" for i in indexes]
    return AtomSelection(
        indexes = indexes, 
        names = names, 
        nmols = nmols, 
        natomspermol = natomspermol,
        groups = groups
    )
end

@testitem "AtomSelection indexes" begin
    using PDBTools
    pdbfile = ComplexMixtures.Testing.pdbfile
    indexes = index.(readPDB(pdbfile, "protein and residue 2"))
    s = AtomSelection(indexes, nmols = 1, natomspermol = 11)
    @test s.imol == ones(Int, 11)
    @test s.index == [12 + i for i = 1:11]
    @test s.names == ["$i" for i in s.index]
    @test s.natoms == 11
    @test s.natomspermol == 11
    @test s.nmols == 1
end

# Function to initialize the structures
function AtomSelection(
    indexes::Vector{Int},
    names::Vector{String};
    nmols::Int = 0,
    natomspermol::Int = 0,
)

    if nmols == 0 && natomspermol == 0
        throw(ArgumentError("Set nmols or natomspermol when defining a selection."))
    end

    natoms = length(indexes)
    if natoms == 0
        throw(ArgumentError("Vector of atom indexes is empty."))
    end
    if length(names) == 0
        throw(ArgumentError(" Vector of atom names was explicitly provided but is empty. "))
    end
    if length(names) != natoms
        throw(ArgumentError("The vector of atom indexes has a different number of elements than the vector of atom names."))
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

    # Setting the vector that contains the index of the molecule of each atom
    imol = Vector{Int}(undef, natoms)
    iat = 0
    for j = 1:nmols
        for i = 1:natomspermol
            iat = iat + 1
            imol[iat] = j
        end
    end

    return AtomSelection(natoms, nmols, natomspermol, indexes, imol, names[1:natomspermol])
end

# AtomSelection show functions
function Base.show(io::IO, s::AtomSelection)
    print(io, "AtomSelection of $(atoms_str(s.natoms)) belonging to $(mol_str(s.nmols)).")
end

@testitem "AtomSelection - argument errors" begin
    using ComplexMixtures
    @test_throws ArgumentError AtomSelection([1,2,3], ["A", "B", "C"])
    @test_throws ArgumentError AtomSelection(Int[], ["A", "B", "C"]; nmols = 1)
    @test_throws ArgumentError AtomSelection([1, 2, 3], String[]; nmols = 1)
    @test_throws ArgumentError AtomSelection([1, 2, 3], String["A", "B"]; nmols = 1)
    @test_throws ArgumentError AtomSelection([1,2,3], ["A", "B", "C"]; nmols = 2)
    @test_throws ArgumentError AtomSelection([1,2,3], ["A", "B", "C"]; natomspermol = 2)
    @test AtomSelection([1,2,3], ["A","B","C"], natomspermol=1) == AtomSelection(3, 3, 1, [1,2,3], [1,2,3], ["A"])
    @test AtomSelection([1,2,3], ["A","B","C"], nmols=1) == AtomSelection(3, 1, 3, [1,2,3], [1,1,1], ["A","B","C"])
end
