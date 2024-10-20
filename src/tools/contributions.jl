"""
    contributions(R::Result, group::Union{SoluteGroup,SolventGroup}; type = :mddf)

Returns the contributions of the atoms of the solute or solvent to the MDDF, coordination number, or MD count.

# Arguments

- `R::Result`: The result of a calculation.
- `group::Union{SoluteGroup,SolventGroup}`: The group of atoms to consider.
- `type::Symbol`: The type of contributions to return. Can be `:mddf` (default), `:coordination_number`, or `:md_count`.

# Examples

```jldoctest
julia> using ComplexMixtures, PDBTools

julia> dir = ComplexMixtures.Testing.data_dir*"/Gromacs";

julia> atoms = readPDB(dir*"/system.pdb");

julia> protein = select(atoms, "protein");

julia> emim = select(atoms, "resname EMI"); 

julia> solute = AtomSelection(protein, nmols = 1)
AtomSelection 
    1231 atoms belonging to 1 molecule(s).
    Atoms per molecule: 1231
    Number of groups: 1231

julia> solvent = AtomSelection(emim, natomspermol = 20)
AtomSelection 
    5080 atoms belonging to 254 molecule(s).
    Atoms per molecule: 20
    Number of groups: 20

julia> results = load(dir*"/protein_EMI.json"); # load pre-calculated results

julia> ca_cb = contributions(results, SoluteGroup(["CA", "CB"])); # contribution of CA and CB atoms to the MDDF

julia> ca_cb = contributions(results, SoluteGroup(["CA", "CB"]); type=:coordination_number); # contribution of CA and CB atoms to the coordination number
```

"""
function contributions(
    R::Result, 
    group::Union{SoluteGroup,SolventGroup};
    type = :mddf, 
    _warn_zero_md_count = true,
    _unsafe_types_from_indices = false,
)

    if !(type in (:mddf, :coordination_number, :md_count))
        throw(ArgumentError("type must be :mddf (default), :coordination_number, or :md_count"))
    end

    if _warn_zero_md_count && type == :mddf && all(==(0), R.md_count_random)
        @warn begin """\n
            This is probably a `coordination_number` calculation only. 
            
            To compute the contributions to coordination numbers, use `type=:coordination_number`.

        """
        end _file=nothing _line=nothing
    end

    if group isa SoluteGroup
        atsel = R.solute
        group_count = R.solute_group_count
    elseif group isa SolventGroup
        atsel = R.solvent
        group_count = R.solvent_group_count
    end

    # If custom groups were provided, we cannot retrieve general group contributions from
    # the group_count array.
    if atsel.custom_groups
        if isnothing(group.group_index) && isnothing(group.group_name)
            throw(ArgumentError("""\n
                Custom groups are defined. Cannot retrieve general group contributions. 
                Please provide a group name or index.
                For example, use SoluteGroup(1) or SoluteGroup("Group1_NAME")
            """)) 
        end
    else
        if isnothing(group.atom_indices) && isnothing(group.atom_names)
            throw(ArgumentError("""\n
                No custom groups are defined. Please provide vectors of *atomic* indices or names.
                For example, to get the contribution of atoms 1, 2 and 3, use SoluteGroup([1,2,3]). 
            """))
        end
    end
    sel_count = zeros(length(group_count[1]))

    # If the index of the groups was provided
    if !isnothing(group.group_index)
        igroup = group.group_index
        if igroup > length(group_count)
            throw(ArgumentError("Group $igroup greater than number of groups ($(length(group_count))) of group contribution array."))
        end
        sel_count .= group_count[igroup]
    end

    # If the name of the group was provided
    if !isnothing(group.group_name)
        igroup = findfirst(==(group.group_name), atsel.group_names)
        if isnothing(igroup) || igroup > length(group_count)
            throw(ArgumentError("Group '$(group.group_name)' not found in the group contribution count array."))
        end
        sel_count .= group_count[igroup]
    end

    # If, instead, atom indices or names were provided, sum over the contributions of the atoms.
    # This sum is different if the structure contanis one or more than one molecule.
    # If the structure has more than one molecule, than the indices are indices of 
    # the atoms *within* the molecule. 

    # Given atom indices, sum over the contributions of the atoms
    if !isnothing(group.atom_indices) 
        if isempty(group.atom_indices)
            throw(ArgumentError("Group selection by group indices is empty"))
        end
        # Check consistency of input indexes
        for i in group.atom_indices
            if all(!=(i), atsel.indices)
                throw(ArgumentError("""\n
                    Atom index $i not found in atom selection. 
                    Indices array: $(print_vector_summary(atsel.indices)).

                    Notes: 
                    - The indices must correspond to the indices of the atoms in the original structure file.
                    - If the atom selection contains more than one molecule, all the atoms corresponding to the 
                      same index *within each molecule* are considered equivalent and summed.

                """
                ))
            end
        end
        # Now run over the types, and sum the contributions. If the selection of 
        # indices have repeated types, the contributions are then *not* summed. 
        if _unsafe_types_from_indices # this only works if: 1) there's only 1 molecule; 2) the indices correspond to the group types
            if atsel.nmols > 1
                throw(ArgumentError("""\n
                    There is more than one molecule in this $(typeof(atsel)). 
                    This is not compatible with `_unsafe_types_from_indices = true`.

                """))
            end
            for i in group.atom_indices
                itype = atom_type(i, atsel.natomspermol; first = atsel.indices[1])
                sel_count .+= group_count[itype]
            end
        else # search for the correspondence between indices and types 
            for itype in eachindex(atsel.group_names)
                if any(iat -> atom_type(iat, atsel.natomspermol; first = atsel.indices[1]) == itype, group.atom_indices)
                    sel_count .+= group_count[itype]
                end
            end
        end
    end

    # Given atom names, sum over the contributions of the atoms
    if !isnothing(group.atom_names) 
        if isempty(group.atom_names)
            throw(ArgumentError("Group selection by group names is empty"))
        end
        # Check consistency of input names
        for name in group.atom_names
            if all(!=(name), atsel.group_names)
                throw(ArgumentError("Atom (or group) name $name not found in group name data."))
            end
        end
        # Now run over the names, and sum the contributions. If the selection of
        # names have repeated types, the contributions are then *not* summed.
        for (itype, name) in enumerate(atsel.group_names)
            if any(==(name), group.atom_names)
                sel_count .+= group_count[itype]
            end
        end
    end

    # Convert to the desired type
    if type == :mddf
        for i in eachindex(R.md_count_random, sel_count)
            if R.md_count_random[i] == 0.0
                sel_count[i] = 0.0
            else
                sel_count[i] /= R.md_count_random[i]
            end
        end
    elseif type == :coordination_number
        sel_count .= cumsum(sel_count)
    elseif type == :md_count
        # do nothing, already md_count
    end

    return sel_count
end

function warning_nmols_types()
    println("""
        WARNING: There is more than one molecule in this selection.
                 Contributions are summed over all atoms of the same type.
    """)
end

@testitem "contributions" begin
    using ComplexMixtures
    using PDBTools: select, readPDB, Select
    using ComplexMixtures.Testing: data_dir
    atoms = readPDB("$data_dir/PDB/trajectory.pdb", "model 1")
    protein = select(atoms, "protein")
    tmao = select(atoms, "resname TMAO")
    solute = AtomSelection(
        protein, nmols = 1;
        group_names = [ "acidic", "basic", "polar", "nonpolar" ],
        group_atom_indices = [ 
            findall(Select("acidic"), protein),
            findall(Select("basic"), protein),
            findall(Select("polar"), protein),
            findall(Select("nonpolar"), protein),
        ]
    )
    solvent = AtomSelection(tmao, natomspermol = 14)
    traj = Trajectory("$data_dir/PDB/trajectory.pdb", solute, solvent, format = "PDBTraj")
    results = mddf(traj)

    @test sum(contributions(results, SoluteGroup("acidic"); type = :md_count)) ≈ 4.4
    @test sum(contributions(results, SoluteGroup("basic"); type = :md_count)) ≈ 2.4
    @test sum(contributions(results, SoluteGroup("polar"); type = :md_count)) ≈ 20.0
    @test sum(contributions(results, SoluteGroup("nonpolar"); type = :md_count)) ≈ 9.4
    
    @test sum(contributions(results, SolventGroup(solvent); type = :md_count)) ≈ 29.4
    @test sum(contributions(results, SolventGroup(tmao); type = :md_count)) ≈ 29.4
    @test sum(contributions(results, SolventGroup(findall(Select("resname TMAO and resnum 1"), atoms)); type = :md_count)) ≈ 29.4

    @test_throws ArgumentError contributions(results, SolventGroup(solvent); type = :wrong_type)
    @test_throws ArgumentError contributions(results, SoluteGroup([1,2,3]))
    @test_throws ArgumentError contributions(results, SoluteGroup(["N", "CA"]))
    @test_throws ArgumentError contributions(results, SolventGroup("acidic"))
    @test_throws ArgumentError contributions(results, SolventGroup(1))
    @test_throws ArgumentError contributions(results, SolventGroup([5000]))
    @test_throws ArgumentError contributions(results, SoluteGroup("NOT_FOUND"))
    @test_throws ArgumentError contributions(results, SoluteGroup(1000))
    @test_throws ArgumentError contributions(results, SolventGroup(Int[]))
    @test_throws ArgumentError contributions(results, SolventGroup(String[]))
    @test_throws ArgumentError contributions(results, SolventGroup(["NOT_FOUND"]))

    solute = AtomSelection(protein, nmols = 1)
    traj = Trajectory("$data_dir/PDB/trajectory.pdb", solute, solvent, format = "PDBTraj")
    results = mddf(traj)
    @test_throws ArgumentError contributions(results, SoluteGroup("acidic"))
    @test_throws ArgumentError contributions(results, SoluteGroup([50000]))

    # Contributions to coordination-number only
    results = coordination_number(traj)

end

@testitem "custom group contributions" begin
    using ComplexMixtures: mddf, Trajectory, Options, AtomSelection, load,
        SoluteGroup, SolventGroup, contributions
    using PDBTools: readPDB, select, Select, iswater, chain, resnum
    using ComplexMixtures.Testing: data_dir

    dir = "$data_dir/NAMD"
    atoms = readPDB("$dir/structure.pdb")
    options = Options(stride = 5, seed = 321, StableRNG = true, nthreads = 1, silent = true)

    #
    # Test computation of custom group counters
    #
    solute = AtomSelection(
        select(atoms , "water and residue 301"), 
        nmols = 1,
        group_atom_indices = [
            findall(Select("water and residue 301 and name H1"), atoms),
            findall(Select("water and residue 301 and name H2"), atoms),
            findall(Select("water and residue 301 and name OH2"), atoms)
        ],
        group_names = ["H1", "H2", "OH2"] 
    )
    solvent = AtomSelection(
        select(atoms, "water and not residue 301"),
        natomspermol = 3, 
        group_atom_indices = [
            findall(Select("water and name H1 and not residue 301"), atoms),
            findall(Select("water and name H2 and not residue 301"), atoms),
            findall(Select("water and name OH2 and not residue 301"), atoms)
        ],
        group_names = ["H1", "H2", "OH2"] 
    )
    traj = Trajectory("$dir/trajectory.dcd", solute, solvent)
    R = mddf(traj)

    @test !all(==(0), contributions(R, SoluteGroup("H1")))
    @test !all(==(0), contributions(R, SoluteGroup("H2")))
    @test !all(==(0), contributions(R, SoluteGroup("OH2")))

    @test !all(==(0), contributions(R, SolventGroup("H1")))
    @test !all(==(0), contributions(R, SolventGroup("H2")))
    @test !all(==(0), contributions(R, SolventGroup("OH2")))

    @test contributions(R, SoluteGroup("H1")) + 
          contributions(R, SoluteGroup("H2")) + 
          contributions(R, SoluteGroup("OH2")) ≈ 
          contributions(R, SolventGroup("H1")) + 
          contributions(R, SolventGroup("H2")) + 
          contributions(R, SolventGroup("OH2"))

    @test contributions(R, SoluteGroup("H1")) + 
          contributions(R, SoluteGroup("H2")) + 
          contributions(R, SoluteGroup("OH2")) ≈ 
          R.mddf

    # Group contributions in autocorrelation computation
    solute = AtomSelection(
        select(atoms, "water and resnum <= 1000"), 
        natomspermol = 3,
        group_atom_indices = [
            findall(Select("water and name H1 and resnum <= 1000"), atoms),
            findall(Select("water and name H2 and resnum <= 1000"), atoms),
            findall(Select("water and name OH2 and resnum <= 1000"), atoms)
        ],
        group_names = ["H1", "H2", "OH2"] 
    )
    traj = Trajectory("$dir/trajectory.dcd", solute)
    R = mddf(traj, options)

    @test !all(==(0), contributions(R, SoluteGroup("H1")))
    @test !all(==(0), contributions(R, SoluteGroup("H2")))
    @test !all(==(0), contributions(R, SoluteGroup("OH2")))

    @test !all(==(0), contributions(R, SolventGroup("H1")))
    @test !all(==(0), contributions(R, SolventGroup("H2")))
    @test !all(==(0), contributions(R, SolventGroup("OH2")))

    @test contributions(R, SoluteGroup("H1")) + 
          contributions(R, SoluteGroup("H2")) + 
          contributions(R, SoluteGroup("OH2")) ≈ 
          contributions(R, SolventGroup("H1")) + 
          contributions(R, SolventGroup("H2")) + 
          contributions(R, SolventGroup("OH2"))

    @test contributions(R, SoluteGroup("H1")) + 
          contributions(R, SoluteGroup("H2")) + 
          contributions(R, SoluteGroup("OH2")) ≈ 
          R.mddf

    # Group contributions in autocorrelation computation, with different chains
    solute = AtomSelection(
        filter(
            at -> iswater(at) && (chain(at) == "B" || chain(at) == "C") && resnum(at) <= 1000, 
            atoms
        ),
        natomspermol = 3, 
        group_atom_indices = [
            findall(Select("water and chain B and resnum <= 1000"), atoms),
            findall(Select("water and chain C and resnum <= 1000"), atoms)
        ],
        group_names = ["B", "C"]
    )
    traj = Trajectory("$dir/trajectory.dcd", solute)
    R = mddf(traj, options)
    @test !all(==(0), contributions(R, SoluteGroup("B")))
    @test !all(==(0), contributions(R, SoluteGroup("C")))
    @test contributions(R, SolventGroup("B")) ≈ contributions(R, SoluteGroup("B"))
    @test contributions(R, SolventGroup("C")) ≈ contributions(R, SoluteGroup("C"))
    @test contributions(R, SoluteGroup("B")) + contributions(R, SoluteGroup("C")) ≈ R.mddf

end

@testitem "shuffled indices custom group contributions" begin
    using ComplexMixtures: mddf, Trajectory, Options, AtomSelection, load,
        SoluteGroup, SolventGroup, contributions
    using PDBTools: readPDB, select, Select, iswater, chain, resnum
    using ComplexMixtures.Testing: data_dir
    import Random: shuffle!

    dir = "$data_dir/NAMD"
    atoms = readPDB("$dir/structure.pdb")
    options = Options(stride = 5, seed = 321, StableRNG = true, nthreads = 1, silent = true)

    #
    # Test computation of custom group counters
    #
    solute = AtomSelection(
        select(atoms , "water and residue 301"), 
        nmols = 1,
        group_atom_indices = [
            shuffle!(findall(Select("water and residue 301 and name H1"), atoms)),
            shuffle!(findall(Select("water and residue 301 and name H2"), atoms)),
            shuffle!(findall(Select("water and residue 301 and name OH2"), atoms))
        ],
        group_names = ["H1", "H2", "OH2"] 
    )
    solvent = AtomSelection(
        select(atoms, "water and not residue 301"),
        natomspermol = 3, 
        group_atom_indices = [
            shuffle!(findall(Select("water and name H1 and not residue 301"), atoms)),
            shuffle!(findall(Select("water and name H2 and not residue 301"), atoms)),
            shuffle!(findall(Select("water and name OH2 and not residue 301"), atoms))
        ],
        group_names = ["H1", "H2", "OH2"] 
    )
    traj = Trajectory("$dir/trajectory.dcd", solute, solvent)
    R = mddf(traj)

    @test !all(==(0), contributions(R, SoluteGroup("H1")))
    @test !all(==(0), contributions(R, SoluteGroup("H2")))
    @test !all(==(0), contributions(R, SoluteGroup("OH2")))

    @test !all(==(0), contributions(R, SolventGroup("H1")))
    @test !all(==(0), contributions(R, SolventGroup("H2")))
    @test !all(==(0), contributions(R, SolventGroup("OH2")))

    @test contributions(R, SoluteGroup("H1")) + 
          contributions(R, SoluteGroup("H2")) + 
          contributions(R, SoluteGroup("OH2")) ≈ 
          contributions(R, SolventGroup("H1")) + 
          contributions(R, SolventGroup("H2")) + 
          contributions(R, SolventGroup("OH2"))

    @test contributions(R, SoluteGroup("H1")) + 
          contributions(R, SoluteGroup("H2")) + 
          contributions(R, SoluteGroup("OH2")) ≈ 
          R.mddf

    # Group contributions in autocorrelation computation
    solute = AtomSelection(
        select(atoms, "water and resnum <= 1000"), 
        natomspermol = 3,
        group_atom_indices = [
            shuffle!(findall(Select("water and name H1 and resnum <= 1000"), atoms)),
            shuffle!(findall(Select("water and name H2 and resnum <= 1000"), atoms)),
            shuffle!(findall(Select("water and name OH2 and resnum <= 1000"), atoms))
        ],
        group_names = ["H1", "H2", "OH2"] 
    )
    traj = Trajectory("$dir/trajectory.dcd", solute)
    R = mddf(traj, options)

    @test !all(==(0), contributions(R, SoluteGroup("H1")))
    @test !all(==(0), contributions(R, SoluteGroup("H2")))
    @test !all(==(0), contributions(R, SoluteGroup("OH2")))

    @test !all(==(0), contributions(R, SolventGroup("H1")))
    @test !all(==(0), contributions(R, SolventGroup("H2")))
    @test !all(==(0), contributions(R, SolventGroup("OH2")))

    @test contributions(R, SoluteGroup("H1")) + 
          contributions(R, SoluteGroup("H2")) + 
          contributions(R, SoluteGroup("OH2")) ≈ 
          contributions(R, SolventGroup("H1")) + 
          contributions(R, SolventGroup("H2")) + 
          contributions(R, SolventGroup("OH2"))

    @test contributions(R, SoluteGroup("H1")) + 
          contributions(R, SoluteGroup("H2")) + 
          contributions(R, SoluteGroup("OH2")) ≈ 
          R.mddf

    # Group contributions in autocorrelation computation, with different chains
    solute = AtomSelection(
        filter(
            at -> iswater(at) && (chain(at) == "B" || chain(at) == "C") && resnum(at) <= 1000, 
            atoms
        ),
        natomspermol = 3, 
        group_atom_indices = [
            shuffle!(findall(Select("water and chain B and resnum <= 1000"), atoms)),
            shuffle!(findall(Select("water and chain C and resnum <= 1000"), atoms))
        ],
        group_names = ["B", "C"]
    )
    traj = Trajectory("$dir/trajectory.dcd", solute)
    R = mddf(traj, options)
    @test !all(==(0), contributions(R, SoluteGroup("B")))
    @test !all(==(0), contributions(R, SoluteGroup("C")))
    @test contributions(R, SolventGroup("B")) ≈ contributions(R, SoluteGroup("B"))
    @test contributions(R, SolventGroup("C")) ≈ contributions(R, SoluteGroup("C"))
    @test contributions(R, SoluteGroup("B")) + contributions(R, SoluteGroup("C")) ≈ R.mddf

end
