#=

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

=#
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

#=
    update_md(md1::MinimumDistance{T}, md2::MinimumDistance{T}) where {T}

$(INTERNAL)

Function that returns the updated minimum distance structure after comparing two structures
associated with the same molecule.

=#
function update_md(md1::MinimumDistance, md2::MinimumDistance)
    ref_atom_within_cutoff = md1.ref_atom_within_cutoff || md2.ref_atom_within_cutoff
    dref = ref_atom_within_cutoff ? min(md1.d_ref_atom, md2.d_ref_atom) : +Inf
    md = if md1.d < md2.d
        MinimumDistance(md1.within_cutoff, md1.i, md1.j, md1.d, ref_atom_within_cutoff, dref)
    else
        MinimumDistance(md2.within_cutoff, md2.i, md2.j, md2.d, ref_atom_within_cutoff, dref)
    end
    return md
end

@testitem "update_md" begin
    import ComplexMixtures as CM
    md1 = CM.MinimumDistance(true, 1, 2, 1.0, true, 1.0)
    md2 = CM.MinimumDistance(true, 1, 2, 0.5, true, 0.5)
    @test CM.update_md(md1, md2) == CM.MinimumDistance(true, 1, 2, 0.5, true, 0.5)
end

#
# Methods to allow multi-threading in CellListMap
#
import CellListMap.PeriodicSystems: copy_output, reset_output!, reducer
copy_output(md::MinimumDistance) = MinimumDistance(
    md.within_cutoff,
    md.i,
    md.j,
    md.d,
    md.ref_atom_within_cutoff,
    md.d_ref_atom,
)
reset_output!(::MinimumDistance) = MinimumDistance(false, 0, 0, +Inf, false, +Inf)
reducer(md1::MinimumDistance, md2::MinimumDistance) = update_md(md1, md2)

#=
    mol_index(i_atom, natomspermol) = (i_atom-1) ÷ natomspermol + 1

# Extended help

Sets the index of the molecule of an atom in the simples situation, in which all 
molecules have the same number of atoms. 

=#
mol_index(i, natomspermol) = (i - 1) ÷ natomspermol + 1

#=
    update_list!(i, j, d2, jref_atom, j_natoms_per_molecule, isolute, list::Vector{MinimumDistance})

Function that updates a list of minimum distances given the indexes of the atoms involved for one pair within cutoff,
for autocorrelations (such that the identity of `isolute` is needed)

=#
function update_list!(
    i,
    j,
    d2,
    jref_atom,
    j_natoms_per_molecule,
    isolute,
    list::Vector{MinimumDistance},
)
    jmol = mol_index(j, j_natoms_per_molecule)
    if jmol != isolute
        d = sqrt(d2)
        ref_atom_within_cutoff = (itype(j, j_natoms_per_molecule) == jref_atom)
        dref = ref_atom_within_cutoff ? d : +Inf
        list[jmol] = update_md(
            list[jmol],
            MinimumDistance(true, i, j, d, ref_atom_within_cutoff, dref),
        )
    end
    return list
end

#=
    update_list!(i, j, d2, jref_atom, j_natoms_per_molecule, list::Vector{MinimumDistance})

Function that updates a list of minimum distances given the indexes of the atoms involved for one pair within cutoff.

=#
function update_list!(
    i,
    j,
    d2,
    jref_atom,
    j_natoms_per_molecule,
    list::Vector{MinimumDistance},
)
    d = sqrt(d2)
    jmol = mol_index(j, j_natoms_per_molecule)
    ref_atom_within_cutoff = (itype(j, j_natoms_per_molecule) == jref_atom)
    dref = ref_atom_within_cutoff ? d : +Inf
    list[jmol] =
        update_md(list[jmol], MinimumDistance(true, i, j, d, ref_atom_within_cutoff, dref))
    return list
end

#=
    minimum_distances!(system::CellListMap.PeriodicSystem, R::Result)

$(INTERNAL)

Function that computes the list of distances of solvent molecules to a solute molecule. 
It updates the lists of minimum distances. 

=#
function minimum_distances!(
    system::AbstractPeriodicSystem,
    R::Result,
    isolute::Int;
    update_lists::Bool,
)
    jref_atom = R.irefatom
    jnatomspermol = R.solvent.natomspermol
    if R.autocorrelation
        map_pairwise!(
            (x, y, i, j, d2, list) ->
                update_list!(i, j, d2, jref_atom, jnatomspermol, isolute, list),
            system;
            update_lists = update_lists,
        )
    else
        map_pairwise!(
            (x, y, i, j, d2, list) ->
                update_list!(i, j, d2, jref_atom, jnatomspermol, list),
            system;
            update_lists = update_lists,
        )
    end
    return system.list
end

#=
    setup_PeriodicSystem(trajectory::Trajectory, options::Options)

$(INTERNAL)

Setup the periodic system from CellListMap, to compute minimimum distances. The system
will be setup such that `xpositions` corresponds to one molecule of the solute, and 
`ypositions` contains all coordinates of all atoms of the solvent. 

=#
function setup_PeriodicSystem(trajectory::Trajectory, options::Options)
    opentraj!(trajectory)
    firstframe!(trajectory)
    nextframe!(trajectory)
    unitcell = convert_unitcell(getunitcell(trajectory)) # returns vector or matrix depending on the data
    closetraj!(trajectory)
    system = PeriodicSystem(
        xpositions = zeros(SVector{3,Float64}, trajectory.solute.natomspermol),
        ypositions = zeros(
            SVector{3,Float64},
            trajectory.solvent.nmols * trajectory.solvent.natomspermol,
        ),
        unitcell = unitcell,
        cutoff = options.usecutoff ? options.cutoff : options.dbulk,
        output = fill(zero(MinimumDistance), trajectory.solvent.nmols),
        output_name = :list,
        lcell = options.lcell,
        parallel = false, # Important: parallellization is performed at the frame level
        autoswap = false, # The lists will be built for the solvent, always
    )
    return system
end

@testitem "setup_PeriodicSystem" begin
    using ComplexMixtures
    using PDBTools
    using ComplexMixtures.Testing
    using StaticArrays
    import CellListMap

    atoms = readPDB(Testing.pdbfile)
    options = Options(stride = 5, seed = 321, StableRNG = true, nthreads = 1, silent = true)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol = 14)

    # Cross-correlation
    protein = AtomSelection(select(atoms, "protein"), nmols = 1)
    traj = Trajectory("$(Testing.data_dir)/NAMD/trajectory.dcd", protein, tmao)
    system = ComplexMixtures.setup_PeriodicSystem(traj, options)
    @test system.cutoff == 10.0
    @test system.list == fill(zero(ComplexMixtures.MinimumDistance), 181)
    @test system.output == fill(zero(ComplexMixtures.MinimumDistance), 181)
    @test system.parallel == false
    @test length(system.xpositions) == 1463
    @test length(system.ypositions) == 2534
    @test system.unitcell ≈ [84.42188262939453 0.0 0.0; 0.0 84.42188262939453 0.0; 0.0 0.0 84.42188262939453]
    @test system._box == CellListMap.Box(
        ComplexMixtures.convert_unitcell(ComplexMixtures.getunitcell(traj)),
        10.0, lcell = options.lcell
    )

    # Auto-correlation
    traj = Trajectory("$(Testing.data_dir)/NAMD/trajectory.dcd", tmao)
    system = ComplexMixtures.setup_PeriodicSystem(traj, options)
    @test system.cutoff == 10.0
    @test system.list == fill(zero(ComplexMixtures.MinimumDistance), 181) # one molecule less
    @test system.output == fill(zero(ComplexMixtures.MinimumDistance), 181)
    @test system.parallel == false
    @test length(system.xpositions) == 14 # one TMAO molecule
    @test length(system.ypositions) == 2534 # one molecule less
    @test system.unitcell ≈ [84.42188262939453 0.0 0.0; 0.0 84.42188262939453 0.0; 0.0 0.0 84.42188262939453]
    @test system._box == CellListMap.Box(
        ComplexMixtures.convert_unitcell(ComplexMixtures.getunitcell(traj)),
        10.0, lcell = options.lcell
    )

end
