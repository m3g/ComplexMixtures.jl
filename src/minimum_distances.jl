"""
    setup_PeriodicSystem(autocorrelation::Bool, trajectory::Trajectory, options::Options)

$(INTERNAL)

Setup the periodic system from CellListMap, to compute minimimum distances. The system
will be setup such that `xpositions` corresponds to one molecule of the solute, and 
`ypositions` contains all coordinates of all atoms of the solvent. This is the same setup
wether an aucorrelation function is being computed, or not.

"""
function setup_PeriodicSystem(autocorrelation::Bool, trajectory::Trajectory, options::Options)
    nextframe!(trajectory)
    unit_cell = setunitcell(trajectory) # returns vector or matrix depending on the data
    firstframe!(trajectory)
    # For autocorrelations, the length of the solvent array contains one less molecule
    ylength = length(trajectory.x_solvent) - R.autocorrelation*trajectory.solvent.n_atoms_per_molecule
    system = CellListMap.PeriodicSystem(
        xpositions = similar(trajectory.x_solute, trajectory.solute.n_atoms_per_molecule),
        ypositions = similar(trajectory.x_solvent, ylength),
        unit_cell = unit_cell,
        cutoff = options.cutoff,
        output = fill(zero(MinimumDistances()), trajectory.solvent.nmols),
        output_name = :list,
        lcell=options.lcell,
        parallel=false, # Important: parallellization is performed at the frame level
    )
    return system
end

"""

$(TYPEDEF)

$(INTERNAL)

# Extended help

This structure contains the information, for each  molecule, of if it is within the 
cutoff distance of the solute, the atom indexes of the associated minimum distance,
the distance, and a label to mark if the reference atom of the molecule is within
the cutoff distance of the solute.

The lists of minimum-distances are stored in arrays of type `Vector{MinimumDistance}`. The index
of this vector corresponds to the index of the molecule in the original array.

$(TYPEDFIELDS)

"""
struct MinimumDistance
    within_cutoff::Bool
    i::Int
    j::Int
    d::Float64
    ref_atom_within_cutoff::Bool
    d_ref_atom::Float64
end
import Base: zero
zero(::Type{MinimumDistance}) = MinimumDistance(false, 0, 0, +Inf, false, +Inf)

"""

```
update_md(md1::MinimumDistance{T}, md2::MinimumDistance{T}) where {T}
```

$(INTERNAL)

Function that returns the updated minimum distance structure after comparing two structures
associated with the same molecule.

"""
function update_md(md1::MinimumDistance, md2::MinimumDistance)
    found_ref = md1.ref_atom_within_cutoff || md2.ref_atom_within_cutoff
    dref = found_ref ? min(md1.d_ref_atom, md2.d_ref_atom) : +Inf
    md = if md1.d < md2.d
        MinimumDistance(md1.within_cutoff, md1.i, md1.j, md1.d, found_ref, dref)
    else
        MinimumDistance(md2.within_cutoff, md2.i, md2.j, md2.d, found_ref, dref)
    end
    return md
end

@testitem "update_md" begin
    import ComplexMixtures as CM
    md1 = cm.MinimumDistance(true, 1, 2, 1.0, true, 1.0)
    md2 = cm.MinimumDistance(true, 1, 2, 0.5, true, 0.5)
    @test CM.update_md(md1,md2) = MinimumDistance(true, 1, 2, 0.5, true, 0.5)
end

#
# Methods to allow multi-threading in CellListMap
#
import CellListMap.PeriodicSystems: copy_output, reset_output!, reducer
copy_output(md::MinimumDistance) = 
    MinimumDistance(md.within_cutoff, md.i, md.j, md.d, md.ref_atom_within_cutoff, md.d_ref_atom)
reset_output!(::MinimumDistance) = 
    MinimumDistance(false, 0, 0, +Inf, false, +Inf)
reducer(md1::MinimumDistance, md2::MinimumDistance) = update_md(md1, md2)

"""
    mol_index(i_atom, n_atoms_per_molecule) = (i_atom-1) ÷ n_atoms_per_molecule + 1

$(INTERNAL)

# Extended help

Sets the index of the molecule of an atom in the simples situation, in which all 
molecules have the same number of atoms. 

"""
mol_index(i, n_atoms_per_molecule) = (i - 1) ÷ n_atoms_per_molecule + 1

"""
    update_list!(i, j, d2, iref_atom::Int, mol_index_i::F, list::Vector{MinimumDistance{T}}) where {F<:Function, T}

$(INTERNAL)

Function that updates a list of minimum distances given the indexes of the atoms involved for one pair within cutoff.

"""
function update_list!(i, j, d2, jref_atom, j_natoms_per_molecule, list::Vector{MinimumDistance})
    d = sqrt(d2)
    jmol = mol_index(j, j_natoms_per_molecule)
    found_ref = j%jref_atom == 0
    dref = found_ref ? d : +Inf
    list[jmol] = update_md(list[jmol], MinimumDistance(true, i, j, d, found_ref, dref))
    return list
end

"""
    minimum_distances!(system::CellListMap.PeriodicSystem)

Function that computes the list of distances of solvent molecules to a solute molecule. 
It updates the lists of minimum distances. 

$(INTERNAL)

# Extended help

"""
function minimum_distances!(system::AbstractPeriodicSystem, options)
    jref_atom = options.irefatom
    jn_atoms_per_molecule = options.solvent.n_atoms_per_molecule
    map_pairwise!(
        (x, y, i, j, d2, list) -> update_list!(i, j, d2, jref_atom, jn_atoms_per_molecule, list),
        system
    )
    return system.list
end
