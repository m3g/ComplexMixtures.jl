"""

$(TYPEDEF)

Structure that contains the information about the solute and solvent molecules.

$(TYPEDFIELDS)

"""
struct Selection
    natoms::Int # Total number of atoms
    nmols::Int # Number of molecules
    natomspermol::Int # Number of atoms per molecule
    index::Vector{Int} # Index of each atom in the full vector of coordinates
    imol::Vector{Int} # index of the molecule to which each atom belongs
    names::Vector{String} # Types of the atoms, to be used in the atom-contributions
end

# Initialize providing the file name, and calling by default PDBTools.select
function Selection(file::String, selection::String; nmols::Int = 0, natomspermol::Int = 0)
    sel = PDBTools.readPDB(file, selection)
    return Selection(sel, nmols = nmols, natomspermol = natomspermol)
end

@testitem "Selection PDBTools" begin
    pdbfile = ComplexMixtures.Testing.pdbfile
    s = Selection(pdbfile, "protein and residue 2", nmols = 1, natomspermol = 11)
    @test s.imol == ones(Int, 11)
    @test s.index == [12 + i for i = 1:11]
    @test s.names == ["N", "HN", "CA", "HA", "CB", "HB1", "HB2", "SG", "HG1", "C", "O"]
    @test s.natoms == 11
    @test s.natomspermol == 11
    @test s.nmols == 1
end

# If the input is a vector of PDBTools.Atom types, load the index and types
function Selection(atoms::Vector{PDBTools.Atom}; nmols::Int = 0, natomspermol::Int = 0)
    indexes = [at.index for at in atoms]
    names = [at.name for at in atoms]
    return Selection(indexes, names, nmols = nmols, natomspermol = natomspermol)
end

@testitem "Selection Vector{PDBTools.Atom}" begin
    import PDBTools
    pdbfile = ComplexMixtures.Testing.pdbfile
    atoms = PDBTools.readPDB(pdbfile, "protein and residue 2")
    s = Selection(atoms, nmols = 1, natomspermol = 11)
    @test s.imol == ones(Int, 11)
    @test s.index == [12 + i for i = 1:11]
    @test s.names == ["N", "HN", "CA", "HA", "CB", "HB1", "HB2", "SG", "HG1", "C", "O"]
    @test s.natoms == 11
    @test s.natomspermol == 11
    @test s.nmols == 1
end

# If no names are provided, just repeat the indexes
function Selection(indexes::Vector{Int}; nmols::Int = 0, natomspermol::Int = 0)
    names = ["$i" for i in indexes]
    return Selection(indexes, names, nmols = nmols, natomspermol = natomspermol)
end

@testitem "Selection indexes" begin
    using PDBTools
    pdbfile = ComplexMixtures.Testing.pdbfile
    indexes = index.(readPDB(pdbfile, "protein and residue 2"))
    s = Selection(indexes, nmols = 1, natomspermol = 11)
    @test s.imol == ones(Int, 11)
    @test s.index == [12 + i for i = 1:11]
    @test s.names == ["$i" for i in s.index]
    @test s.natoms == 11
    @test s.natomspermol == 11
    @test s.nmols == 1
end

# Function to initialize the structures
function Selection(
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

    return Selection(natoms, nmols, natomspermol, indexes, imol, names[1:natomspermol])
end

# Selection of show functions
function Base.show(io::IO, s::Selection)
    print(io, "Selection of $(atoms_str(s.natoms)) belonging to $(mol_str(s.nmols)).")
end
