"""
    mol_range(imol, n_atoms_per_molecule)

$(INTERNAL)

Given the index and the number of atoms per molecule, returns the range of indices 
of of an array of coordinates that corresponds to the molecule.

"""
function mol_range(imol, n_atoms_per_molecule)
    first = (imol - 1) * n_atoms_per_molecule + 1
    last = first + n_atoms_per_molecule - 1
    return first:last
end

"""
    viewmol(i::Int, x::Vector{T}, n::Int) where T

$(INTERNAL)

Returns a view of a coordinate vector corresponding to the atoms of a molecule with index i. n is the number of atoms of the molecule.

"""
viewmol(i::Int, x::Vector{T}, n::Int) where {T} = @view(x[mol_range(i, n)])

# From the selection of the solute or solvent
viewmol(i::Int, x::Vector{T}, s::Union{Selection,SolSummary}) where {T} =
    viewmol(i, x, s.natomspermol)
