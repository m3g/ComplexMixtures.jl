"""

$(TYPEDEF)

$(INTERNAL)

# Extended help

This structure contains the information, for each  molecule, of if it is within the 
cutoff distance of the solute, the atom indexes of the associated minimum distance,
the distance, and a label to mark if the reference atom of the molecule is within
the cutoff distance of the solute.

The lists of minimum-distances are stored in arrays of type `Vector{MinimumDistance{T}}`. The index
of this vector corresponds to the index of the molecule in the original array.

$(TYPEDFIELDS)

"""
struct MinimumDistance{T}
    within_cutoff::Bool
    i::Int
    j::Int
    d::T
    ref_atom_within_cutoff::Bool
    d_ref_atom::T
end
import Base: zero
zero(::Type{MinimumDistance{T}}) where {T} = MinimumDistance(false, 0, 0, typemax(T), false, typemax(T))

"""

```
update_md(md1::MinimumDistance{T}, md2::MinimumDistance{T}) where {T}
```

$(INTERNAL)

Function that returns the updated minimum distance structure after comparing two structures
associated with the same molecule.


"""
function update_md(md1::MinimumDistance{T}, md2::MinimumDistance{T}) where {T}
    found_ref = md1.ref_atom_within_cutoff || md2.ref_atom_within_cutoff
    dref = found_ref ? min(md1.d_ref_atom, md2.d_ref_atom) : typemax(T)
    if md1.d < md2.d
        md = MinimumDistance{T}(md1.within_cutoff, md1.i, md1.j, md1.d, found_ref, dref)
    else
        md = MinimumDistance{T}(md2.within_cutoff, md2.i, md2.j, md2.d, found_ref, dref)
    end
    return md
end

#
# Methods to allow multi-threading in CellListMap
#
import CellListMap.PeriodicSystems: copy_output, reset_output!, reducer
copy_output(md::MinimumDistance{T}) where {T} = 
    MinimumDistance{T}(md.within_cutoff, md.i, md.j, md.d, md.ref_atom_within_cutoff, md.d_ref_atom)
reset_output!(md1::MinimumDistance{T}) where {T} = 
    MinimumDistance{T}(false, 0, 0, typemax(T), false, typemax(T))
reducer(md1::MinimumDistance{T}, md2::MinimumDistance{T}) where {T} = update_md(md1, md2)

"""

```
mol_index(i_atom,n_atoms_per_molecule) = (i_atom-1) ÷ n_atoms_per_molecule + 1
```

$(INTERNAL)

# Extended help

Sets the index of the molecule of an atom in the simples situation, in which all 
molecules have the same number of atoms. 

"""
mol_index(i, n_atoms_per_molecule) = (i - 1) ÷ n_atoms_per_molecule + 1

"""

$(INTERNAL)

"""
function update_list!(i, j, d2, iref_atom::Int, mol_index_i::F, list::Vector{MinimumDistance{T}},) where {F<:Function, T}
    d = sqrt(d2)
    imol = mol_index_i(i)
    found_ref = i%iref_atom == 0
    dref = found_ref ? d : typemax(T)
    list[imol] = update_md(list[imol], MinimumDistance{T}(true, i, j, d, found_ref, dref))
    return list
end

"""

```
minimum_distances!
```

Function that computes the list of distances of solvent molecules to a solute molecule. 
It updates the lists of minimum distances. 

$(INTERNAL)

# Extended help

```
function minimum_distances!(
    mol_index_i::F,
    x_list::AbstractVector{<:MinimumDistance},
    y::AbstractVector,
    box, cl;
    parallel = true
) where {F<:Function}
```

Compute the list of minimum distances given the precomputed cell lists, auxiliary vectors, and 
the functions that associate each atom to the corresponding molecule. Should be preferred for multiple calls of this function,
with the outer update of the cell lists.

## Example

```
julia> using CellListMap, StaticArrays

julia> x = [ rand(SVector{3,Float64}) for _ in 1:12 ];

julia> y = [ rand(SVector{2,Float64}) for _ in 1:800 ];

julia> box = Box([1,1],0.2);

julia> cl = CellList(x,y,box);

julia> list = fill(zero(MinimumDistance{Float64}), 3) # 3 molecules
3-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)
 MinimumDistance{Float64}(false, 0, 0, Inf)

julia> list_threaded = [ copy(list) for _ in 1:nbatches(cl) ];

julia> minimum_distances!(i -> mol_index(i,4), 1, list, box, cl)
4-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 1, 509, 0.011173109803708009, true, 0.011173109803708009)
 MinimumDistance{Float64}(true, 5, 494, 0.026376853825789335, true, 0.026376853825789335)
 MinimumDistance{Float64}(true, 8, 245, 0.009718966033613377, true, 0.009718966033613377)

```

"""
function minimum_distances!(
    mol_index_i::F,
    iref_atom::Int, 
    x_list::Vector{<:MinimumDistance},
    box::Box,
    cl::CellListMap.CellListPair;
    parallel = true,
    list_threaded = nothing, 
) where {F<:Function}
    reset!(x_list)
    reset!(list_threaded)
    map_pairwise!(
        (x, y, i, j, d2, x_list) -> update_list!(i, j, d2, iref_atom, mol_index_i, x_list),
        x_list,
        box,
        cl;
        parallel = parallel,
        reduce = reduce_list!,
        output_threaded = list_threaded
    )
    return x_list
end

#
# Functions used to reduce the lists in parallel calculations
#
"""

```
reduce_list!(list, list_threaded)
```

$(INTERNAL)

Update the final `list` of minimum-distances given the threaded list `list_threaded`.

"""
function reduce_list!(list, list_threaded)
    @. list = list_threaded[1]
    for it in 2:length(list_threaded)
        @. list = update_md(list, list_threaded[it])
    end
    return list
end

#
# Reset list functions
#
reset!(::Nothing) = nothing
reset!(list::Vector{MinimumDistance{T}}) where {T} = fill!(list, zero(MinimumDistance{T}))
function reset!(list_threaded::Vector{Vector{MinimumDistance{T}}}) where {T}
    for i in eachindex(list_threaded)
        reset!(list_threaded[i])
    end
    return nothing
end

