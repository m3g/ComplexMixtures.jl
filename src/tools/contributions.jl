#
# Error functions for wrong input of atom indices or names
#
function _index_not_found_error(i, atsel)
    throw(ArgumentError("""\n
            Atom index $i not found in atom selection. 
            Indices array: $(print_vector_summary(atsel.indices)).

            Notes: 
            - The indices must correspond to the indices of the atoms in the original structure file.
            - If the atom selection contains more than one molecule, all the atoms corresponding to the 
              same index *within each molecule* are considered equivalent and summed.

        """))
end
function _name_not_found_error(name, atsel)
    throw(ArgumentError("""\n
            Atom (or group) name $name not found in atom selection. 
            Group names list: $(print_vector_summary(atsel.group_names)).

        """))
end


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

julia> dir = ComplexMixtures.data_dir*"/Gromacs";

julia> atoms = read_pdb(dir*"/system.pdb");

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
    type=:mddf,
    _warn_zero_md_count=true,
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
            throw(ArgumentError("""\n
                Group $igroup greater than number of groups ($(length(group_count))) of group contribution array.

            """))
        end
        sel_count .= group_count[igroup]
    end

    # If the name of the group was provided
    if !isnothing(group.group_name)
        igroup = findfirst(==(group.group_name), atsel.group_names)
        isnothing(igroup) && _name_not_found_error(group.group_name, atsel)
        sel_count .= group_count[igroup]
    end

    # If, instead, atom indices or names were provided, sum over the contributions of the atoms.
    # This sum is different if the structure contanis one or more than one molecule.
    # If the structure has more than one molecule, than the indices are indices of 
    # the atoms *within* the molecule. 

    # Given atom indices, sum over the contributions of the atoms
    if !isnothing(group.atom_indices)
        if isempty(group.atom_indices)
            throw(ArgumentError("""\n 
                Group selection by group indices is empty.

            """))
        end
        if !allunique(group.atom_indices)
            throw(ArgumentError("""\n
                Selection by atom indices contains repeated indices.

            """))
        end
        # For a single molecule, all contributions are summed up
        if atsel.nmols == 1
            for iat in group.atom_indices
                itype = findfirst(==(iat), atsel.indices)
                isnothing(itype) && _index_not_found_error(iat, atsel)
                sel_count .+= group_count[itype]
            end
        else
            # If there's more than one molecule, the contributions are stored by 
            # atom type. Thus, if the user asked for the contribution of *one* molecule
            # the contributions are scaled accordingly. 
            for iat in group.atom_indices
                any(==(iat), atsel.indices) || _index_not_found_error(iat, atsel)
                itype = atom_type(iat, atsel.natomspermol; first=first(atsel.indices))
                sel_count .+= group_count[itype] / atsel.nmols
            end
        end
    end

    # Given atom names, sum over the contributions of the atoms
    if !isnothing(group.atom_names)
        if isempty(group.atom_names)
            throw(ArgumentError("""\n
                Selection by atom names is empty.

            """))
        end
        if !allunique(group.atom_names)
            throw(ArgumentError("""\n
                Selection by atom names contains repeated names.

            """))
        end
        for atom_name in group.atom_names
            found_atom_name = false
            for (igroup, name) in enumerate(atsel.group_names)
                if atom_name == name
                    found_atom_name = true
                    sel_count .+= group_count[igroup]
                end
            end
            found_atom_name || _name_not_found_error(atom_name, atsel)
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

@testitem "contributions" begin
    using ComplexMixtures
    using ComplexMixtures: data_dir
    using PDBTools
    atoms = read_pdb("$data_dir/PDB/trajectory.pdb", "model 1")
    protein = select(atoms, "protein")
    tmao = select(atoms, "resname TMAO")
    solute = AtomSelection(
        protein, nmols=1;
        group_names=["acidic", "basic", "polar", "nonpolar"],
        group_atom_indices=[
            findall(Select("acidic"), protein),
            findall(Select("basic"), protein),
            findall(Select("polar and not resname GLY and not resname CYS"), protein),
            findall(Select("nonpolar or resname GLY or resname CYS"), protein),
        ]
    )
    solvent = AtomSelection(tmao, natomspermol=14)
    traj = Trajectory("$data_dir/PDB/trajectory.pdb", solute, solvent, format="PDBTraj")
    results = mddf(traj)

    @test sum(contributions(results, SoluteGroup("acidic"); type=:md_count)) ≈ 4.4
    @test sum(contributions(results, SoluteGroup("basic"); type=:md_count)) ≈ 2.4
    @test sum(contributions(results, SoluteGroup("polar"); type=:md_count)) ≈ 20.0
    @test sum(contributions(results, SoluteGroup("nonpolar"); type=:md_count)) ≈ 9.4
    @test sum(contributions(results, SolventGroup(findall(Select("resname TMAO and resnum 1"), atoms)); type=:md_count)) ≈ 29.4 / length(eachresidue(tmao))

    @test_throws ArgumentError contributions(results, SolventGroup(solvent); type=:wrong_type)
    @test_throws ArgumentError contributions(results, SoluteGroup([1, 2, 3]))
    @test_throws ArgumentError contributions(results, SoluteGroup(["N", "CA"]))
    @test_throws ArgumentError contributions(results, SoluteGroup(["N", "N"]))
    @test_throws ArgumentError contributions(results, SolventGroup("acidic"))
    @test_throws ArgumentError contributions(results, SolventGroup(1))
    @test_throws ArgumentError contributions(results, SolventGroup([5000]))
    @test_throws ArgumentError contributions(results, SoluteGroup("NOT_FOUND"))
    @test_throws ArgumentError contributions(results, SoluteGroup(1000))
    @test_throws ArgumentError contributions(results, SolventGroup(Int[]))
    @test_throws ArgumentError contributions(results, SolventGroup(String[]))
    @test_throws ArgumentError contributions(results, SolventGroup(["NOT_FOUND"]))

    solute = AtomSelection(protein, nmols=1)
    traj = Trajectory("$data_dir/PDB/trajectory.pdb", solute, solvent, format="PDBTraj")
    results = mddf(traj)
    @test_throws ArgumentError contributions(results, SoluteGroup("acidic"))
    @test_throws ArgumentError contributions(results, SoluteGroup([50000]))
    @test sum(contributions(results, SolventGroup([1483]))) ≈ 0.3961968338913652
    @test_throws ArgumentError contributions(results, SolventGroup([1483, 1483]))

    # Testing decomposition of solute contributions
    @test results.mddf ≈ contributions(results, SoluteGroup(protein))
    r = collect(eachresidue(protein))
    @test results.mddf ≈ mapreduce(x -> contributions(results, SoluteGroup(x)), +, r)

    # Now test if the solute is a dicontinuous set of atoms in the original structure
    # The solute has only one molecule.
    ala_residues = select(atoms, "resname ALA")
    solute = AtomSelection(ala_residues, nmols=1)
    traj = Trajectory("$data_dir/PDB/trajectory.pdb", solute, solvent, format="PDBTraj")
    results = mddf(traj)
    @test results.mddf ≈ contributions(results, SoluteGroup(ala_residues))
    r = collect(eachresidue(ala_residues))
    @test results.mddf ≈ mapreduce(x -> contributions(results, SoluteGroup(x)), +, r)

    # Contributions of all solvent residues: with multiple molecules, the contributions should
    # be summed over all atoms of the same type, not repeating types. 
    @test results.mddf ≈ contributions(results, SolventGroup(tmao))
    r = collect(eachresidue(tmao))
    @test all(results.mddf ≈ length(r) * contributions(results, SolventGroup(x)) for x in r)

end

@testitem "custom group contributions" begin
    using ComplexMixtures: mddf, Trajectory, Options, AtomSelection, load,
        SoluteGroup, SolventGroup, contributions
    using ComplexMixtures: data_dir
    using PDBTools: read_pdb, select, Select, iswater, chain, resnum

    dir = "$data_dir/NAMD"
    atoms = read_pdb("$dir/structure.pdb")
    options = Options(stride=5, seed=321, StableRNG=true, nthreads=1, silent=true)

    #
    # Test computation of custom group counters
    #
    solute = AtomSelection(
        select(atoms, "water and residue 301"),
        nmols=1,
        group_atom_indices=[
            findall(Select("water and residue 301 and name H1"), atoms),
            findall(Select("water and residue 301 and name H2"), atoms),
            findall(Select("water and residue 301 and name OH2"), atoms)
        ],
        group_names=["H1", "H2", "OH2"]
    )
    solvent = AtomSelection(
        select(atoms, "water and not residue 301"),
        natomspermol=3,
        group_atom_indices=[
            findall(Select("water and name H1 and not residue 301"), atoms),
            findall(Select("water and name H2 and not residue 301"), atoms),
            findall(Select("water and name OH2 and not residue 301"), atoms)
        ],
        group_names=["H1", "H2", "OH2"]
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
        natomspermol=3,
        group_atom_indices=[
            findall(Select("water and name H1 and resnum <= 1000"), atoms),
            findall(Select("water and name H2 and resnum <= 1000"), atoms),
            findall(Select("water and name OH2 and resnum <= 1000"), atoms)
        ],
        group_names=["H1", "H2", "OH2"]
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
        natomspermol=3,
        group_atom_indices=[
            findall(Select("water and chain B and resnum <= 1000"), atoms),
            findall(Select("water and chain C and resnum <= 1000"), atoms)
        ],
        group_names=["B", "C"]
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
    using ComplexMixtures: data_dir
    using PDBTools: read_pdb, select, Select, iswater, chain, resnum
    import Random: shuffle!

    dir = "$data_dir/NAMD"
    atoms = read_pdb("$dir/structure.pdb")
    options = Options(stride=5, seed=321, StableRNG=true, nthreads=1, silent=true)

    #
    # Test computation of custom group counters
    #
    solute = AtomSelection(
        select(atoms, "water and residue 301"),
        nmols=1,
        group_atom_indices=[
            shuffle!(findall(Select("water and residue 301 and name H1"), atoms)),
            shuffle!(findall(Select("water and residue 301 and name H2"), atoms)),
            shuffle!(findall(Select("water and residue 301 and name OH2"), atoms))
        ],
        group_names=["H1", "H2", "OH2"]
    )
    solvent = AtomSelection(
        select(atoms, "water and not residue 301"),
        natomspermol=3,
        group_atom_indices=[
            shuffle!(findall(Select("water and name H1 and not residue 301"), atoms)),
            shuffle!(findall(Select("water and name H2 and not residue 301"), atoms)),
            shuffle!(findall(Select("water and name OH2 and not residue 301"), atoms))
        ],
        group_names=["H1", "H2", "OH2"]
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
        natomspermol=3,
        group_atom_indices=[
            shuffle!(findall(Select("water and name H1 and resnum <= 1000"), atoms)),
            shuffle!(findall(Select("water and name H2 and resnum <= 1000"), atoms)),
            shuffle!(findall(Select("water and name OH2 and resnum <= 1000"), atoms))
        ],
        group_names=["H1", "H2", "OH2"]
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
        natomspermol=3,
        group_atom_indices=[
            shuffle!(findall(Select("water and chain B and resnum <= 1000"), atoms)),
            shuffle!(findall(Select("water and chain C and resnum <= 1000"), atoms))
        ],
        group_names=["B", "C"]
    )
    traj = Trajectory("$dir/trajectory.dcd", solute)
    R = mddf(traj, options)
    @test !all(==(0), contributions(R, SoluteGroup("B")))
    @test !all(==(0), contributions(R, SoluteGroup("C")))
    @test contributions(R, SolventGroup("B")) ≈ contributions(R, SoluteGroup("B"))
    @test contributions(R, SolventGroup("C")) ≈ contributions(R, SoluteGroup("C"))
    @test contributions(R, SoluteGroup("B")) + contributions(R, SoluteGroup("C")) ≈ R.mddf

end
