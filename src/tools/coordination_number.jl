"""
    coordination_number(R::Result) = R.coordination_number
    coordination_number(R::Result, s::Union{SoluteGroup,SolventGroup})

Here, the coordination number is the number of molecules of the solvent that are within a certain distance from the solute.

If no group is defined (first call above), the coordination number of the whole solute or solvent is returned.

If a group is defined (second call above), the contributions of this group of atoms to the coordination number are 
returned, as a function of the distance to the solute: 

- `R` are the results obtained, that is, a `Result` data structure,
- `s` is the solute or solvent selection (of type `ComplexMixtures.AtomSelection`)

# Examples

## Coordination number of a subgroup of atoms of the solute

Here we compute the coordination number of the atoms of the Alanine residues of a protein,
relative to the solvent (TMAO), as a function of the distance. For simplicity, the 
coordination number at 5 Å is shown. We use, in this case, the `SoluteGroup` type to define the group of atoms,
providing `select(ats, "resname ALA")` as the selection of atoms of the solute. `select` is a 
function from the `PDBTools` package. 

```jldoctest
julia> using ComplexMixtures, PDBTools

julia> using ComplexMixtures: data_dir

julia> ats = read_pdb(joinpath(data_dir,"NAMD/structure.pdb"));

julia> solute = AtomSelection(select(ats, "protein"); nmols=1);

julia> solvent = AtomSelection(select(ats, "resname TMAO"); natomspermol=14);

julia> R = mddf(
           joinpath(data_dir,"NAMD/trajectory.dcd"), 
           solute, 
           solvent, 
           Options(bulk_range=(8.0,12.0), silent=true)
        );

julia> i5 = findfirst(>=(5), R.d) # index for distance ≈ 5 Å
251

julia> ala_residues = select(ats, "resname ALA");

julia> coordination_number(R, SoluteGroup(ala_residues))[i5]
0.45
```

Alternatively to the use of the `PDBTools.select` function, an array with the indices of the
atoms of the protein can be used to define the solute group, as, for example:

```julia-repl
julia> ala_indices = findall(at -> resname(at) == "ALA", ats);

julia> coordination_number(R, SoluteGroup(ala_indices))[i5]
0.45
```

A similar syntax can be used to compute contributions of the solvent atoms to the MDDF.

## Coordination numbers if the groups are predefined

If group contributions were precomputed, the name of the group can be used to compute the coordination number:

```jldoctest
julia> using ComplexMixtures, PDBTools

julia> using ComplexMixtures: data_dir

julia> ats = read_pdb(joinpath(data_dir,"NAMD/structure.pdb"));

julia> solute = AtomSelection(
           select(ats, "protein"); 
           nmols=1,
           # predefinition of the group of atoms:
           group_atom_indices = [ findall(at -> resname(at) == "ALA", ats) ],
           group_names = ["alanine residues"]
        )
AtomSelection 
    1463 atoms belonging to 1 molecule(s).
    Atoms per molecule: 1463
    Number of groups: 1

julia> solvent = AtomSelection(select(ats, "resname TMAO"); natomspermol=14);

julia> R = mddf(
           joinpath(data_dir,"NAMD/trajectory.dcd"), 
           solute, 
           solvent, 
           Options(bulk_range=(8.0,12.0), silent=true)
        );

julia> i5 = findfirst(>=(5), R.d) # index for distance ≈ 5 Å
251

julia> coordination_number(R, SoluteGroup("alanine residues"))[i5]
0.45
```

"""
function coordination_number end
coordination_number(R::Result) = R.coordination_number
coordination_number(R::Result, atsel::Union{SoluteGroup,SolventGroup}) = contributions(R, atsel; type=:coordination_number)

@testitem "coordination_number" begin
    using ComplexMixtures: coordination_number, contributions, mddf, Trajectory, Options, AtomSelection, load
    using PDBTools: read_pdb, select
    using ComplexMixtures: data_dir

    dir = "$data_dir/NAMD"
    atoms = read_pdb("$dir/structure.pdb")
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    options = Options(lastframe=1, seed=321, StableRNG=true, nthreads=1, silent=true, n_random_samples=200)
    traj = Trajectory("$dir/trajectory.dcd", protein, tmao)
    R = mddf(traj, options)

    # Test self-consistency
    @test sum(sum.(R.solute_group_count)) ≈ sum(sum.(R.solvent_group_count))
    @test coordination_number(R) == R.coordination_number

    # Checked with vmd: same residue as (resname TMAO and within 3.0 of protein)
    @test R.coordination_number[findfirst(>(3), R.d)] == 7.0
    @test R.coordination_number[findfirst(>(5), R.d)] == 14.0

    # THR93 forms a hydrogen bond with TMAO in this frame
    thr93 = select(atoms, "protein and residue 93")
    @test last(coordination_number(R, SoluteGroup(thr93))) == 1

    # Test coordination number of a solvent atom
    @test sum(coordination_number(R, SolventGroup(["O1"]))) == 1171.0
end
