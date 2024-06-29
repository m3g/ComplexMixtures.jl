"""
    ResidueContributions (data structure)

## Constructor function: 

    ResidueContributions(
        results::Result, atoms::AbstractVector{PDBTools.Atom};
        dmin=1.5, dmax=3.5,
        type=:mddf,
    )

Compute the residue contributions to the solute-solvent pair distribution function.
The function returns a `ResidueContributions` object that can be used to plot the residue contributions,
or to perform arithmetic operations with other `ResidueContributions` objects.

# Arguments

- `results::Result`: The result object obtained from the `mddf` function.
- `atoms::AbstractVector{PDBTools.Atom}`: The vector of atoms of the solute, or a part of it.

# Optional arguments

- `dmin::Float64`: The minimum distance to consider. Default is `1.5`.
- `dmax::Float64`: The maximum distance to consider. Default is `3.5`.
- `type::Symbol`: The type of the pair distribution function (`:mddf`, `:md_count`, or `:coordination_number`). Default is `:mddf`.

A structure of type `ResultContributions` can be used to plot the residue contributions to the solute-solvent pair distribution function,
using the `Plots.contourf` function, and to perform arithmetic operations with other `ResidueContributions` objects.

# Examples

## Constructing a ResidueContributions object

```julia-repl
julia> using ComplexMixtures, PDBTools

julia> using ComplexMixtures.Testing

julia> atoms = readPDB(pdbfile);

julia> protein = AtomSelection(select(atoms, "protein"), nmols=1);

julia> water = AtomSelection(select(atoms, "water"), natomspermol=3);

julia> traj = Trajectory("\$data_dir/NAMD/trajectory.dcd", protein, water);

julia> results = mddf(traj, Options(bulk_range=(8.0, 12.0)))

julia> rc = ResidueContributions(result, select(atoms, "protein"))
         Residue Contributions
    3.51 ████████████████████████████████████████████████████
    3.27 ███████████   ████████ ████████████████ ████████████
    3.03 █████ █████   ████████    █████████████ ████████████
    2.79 █████ █████   ████████    █████████████ ████████████
 d  2.55 █████ █████   ████████    █████████████ ████████████
    2.31 █████ █████   ████████    █████████████ ████████████
    2.07 █████ █████   ████████    █████████████ ████████████
    1.83 ███████████████████████████████████████ ████████████
    1.59 ████████████████████████████████████████████████████
         A1      S17     V33     D49     G65     N81     G97     
```

## Plotting 

```julia
using ComplexMixtures, PDBTools, Plots
...
result = mddf(traj, options)
rc = ResidueContributions(result, select(atoms, "protein"))
contourf(rc) # plots a contour map
```

## Arithmetic operations

```julia
using ComplexMixtures, PDBTools, Plots
...
# first simulation (for example, low temperature):
result1 = mddf(traj2, options)
rc1 = ResidueContributions(result1, select(atoms, "protein"))
# second simulation (for example, high temperature):
result2 = mddf(traj2, options)
rc2 = ResidueContributions(result2, select(atoms, "protein"))
# difference of the residue contributions between the two simulations:
rc_diff = rc2 - rc1
contourf(rc_diff) # plots a contour map of the difference
```

"""
struct ResidueContributions
    d::Vector{Float64}
    xticks::Tuple{Vector{Int},Vector{String}}
    residue_contributions::Matrix{Float64}
end

function ResidueContributions(
    results::Result, atoms::AbstractVector{PDBTools.Atom};
    dmin=1.5, dmax=3.5,
    type=:mddf,
)

    if length(atoms) == 0
        throw(ArgumentError("""\n
            The atom selection provided is empty.
        """))
    end
    if results.solute.custom_groups
        _custom_group_error_for_ResidueContributions()
    end

    # collect the list of residues (using PDBTools)
    residues = collect(PDBTools.eachresidue(atoms))

    # Create matrix that will cotain the contribution per 
    # residue as a function of the distance:
    # number of rows of the mddf histogram is (length(results.d)) and 
    # number of columns equal to the number of residues
    rescontrib = zeros(length(results.d), length(residues))

    # Each column is then filled up with the contributions of each residue
    for (ires, residue) in enumerate(residues)
        rescontrib[:, ires] .= contributions(results, SoluteGroup(residue); type)
    end

    # Plot only for distances within 1.5 and 3.5:
    # find the indexes of the distances that are within the range
    idmin = findfirst(d -> d > dmin, results.d)
    idmax = findfirst(d -> d > dmax, results.d)
    isnothing(idmin) && (idmin = 1)
    isnothing(idmax) && (idmax = length(results.d))

    # Obtain pretty labels for the residues in the x-axis (using PDBTools)
    xticks = PDBTools.residue_ticks(atoms;
        first=PDBTools.resnum(first(residues)),
        last=PDBTools.resnum(last(residues)),
        oneletter=false,
        serial=false,
    )

    return ResidueContributions(
        results.d[idmin:idmax],
        xticks,
        rescontrib[idmin:idmax, 1:length(residues)]
    )
end

function _check_identity_of_residues(rc1::ResidueContributions, rc2::ResidueContributions)
    if rc1.xticks != rc2.xticks
        throw(ArgumentError("The residues in the two ResidueContributions objects differ."))
    end
    if rc1.d != rc2.d
        throw(ArgumentError("The distances in the two ResidueContributions objects differ."))
    end
    return nothing
end

import Base: -, +, /, *
function -(rc1::ResidueContributions, rc2::ResidueContributions)
    _check_identity_of_residues(rc1, rc2)
    return ResidueContributions(rc1.d, rc1.xticks, rc1.residue_contributions - rc2.residue_contributions)
end

function +(rc1::ResidueContributions, rc2::ResidueContributions)
    _check_identity_of_residues(rc1, rc2)
    return ResidueContributions(rc1.d, rc1.xticks, rc1.residue_contributions + rc2.residue_contributions)
end

function /(rc1::ResidueContributions, rc2::ResidueContributions)
    _check_identity_of_residues(rc1, rc2)
    if any(==(0), rc2.residue_contributions)
        @warn begin
            """\n
        Division by zero detected. These elements will be set to zero."

    """
        end _file = nothing _line = nothing
    end
    return ResidueContributions(
        rc1.d,
        rc1.xticks,
        @. ifelse(rc2.residue_contributions == 0, 0.0, rc1.residue_contributions ./ rc2.residue_contributions)
    )

end

function *(rc1::ResidueContributions, rc2::ResidueContributions)
    _check_identity_of_residues(rc1, rc2)
    return ResidueContributions(rc1.d, rc1.xticks, rc1.residue_contributions .* rc2.residue_contributions)
end

const tempo_color_scheme = [
    "#FFF6F4", "#FDF5F3", "#FCF4F1", "#FBF3F0", "#F9F2EE", "#F8F1ED", "#F7F0EB", "#F5EFEA", "#F4EEE8", "#F2EDE7", "#F1ECE5", "#F0EBE4",
    "#EEEAE2", "#EDEAE1", "#EBE9DF", "#EAE8DE", "#E9E7DD", "#E7E6DB", "#E6E5DA", "#E4E4D8", "#E3E3D7", "#E2E2D6", "#E0E2D4", "#DFE1D3",
    "#DDE0D1", "#DCDFD0", "#DBDECF", "#D9DDCD", "#D8DDCC", "#D6DCCB", "#D5DBC9", "#D3DAC8", "#D2D9C7", "#D1D8C5", "#CFD8C4", "#CED7C3",
    "#CCD6C1", "#CBD5C0", "#C9D4BF", "#C8D4BE", "#C6D3BC", "#C5D2BB", "#C3D1BA", "#C2D1B9", "#C0D0B7", "#BFCFB6", "#BDCEB5", "#BCCEB4",
    "#BACDB3", "#B9CCB2", "#B7CBB0", "#B6CBAF", "#B4CAAE", "#B3C9AD", "#B1C8AC", "#B0C8AB", "#AEC7AA", "#ACC6A9", "#ABC5A8", "#A9C5A6",
    "#A8C4A5", "#A6C3A4", "#A4C3A3", "#A3C2A2", "#A1C1A1", "#A0C0A0", "#9EC09F", "#9CBF9F", "#9BBE9E", "#99BE9D", "#97BD9C", "#96BC9B",
    "#94BC9A", "#92BB99", "#91BA98", "#8FBA97", "#8DB997", "#8BB896", "#8AB795", "#88B794", "#86B693", "#85B593", "#83B592", "#81B491",
    "#7FB390", "#7DB390", "#7CB28F", "#7AB18E", "#78B18E", "#76B08D", "#74AF8D", "#72AF8C", "#71AE8B", "#6FAD8B", "#6DAD8A", "#6BAC8A",
    "#69AB89", "#67AB89", "#65AA88", "#63A988", "#61A987", "#5FA887", "#5DA786", "#5BA686", "#59A685", "#57A585", "#56A485", "#54A484",
    "#52A384", "#50A284", "#4EA183", "#4BA183", "#49A083", "#479F82", "#459F82", "#439E82", "#419D82", "#3F9C81", "#3D9C81", "#3B9B81",
    "#3A9A81", "#389981", "#369880", "#349880", "#329780", "#309680", "#2E9580", "#2C947F", "#2A937F", "#29937F", "#27927F", "#25917F",
    "#24907F", "#228F7E", "#218E7E", "#1F8D7E", "#1E8D7E", "#1C8C7E", "#1B8B7D", "#1A8A7D", "#19897D", "#17887D", "#16877C", "#16867C",
    "#15857C", "#14847C", "#13847B", "#13837B", "#12827B", "#12817B", "#11807A", "#117F7A", "#117E7A", "#117D79", "#117C79", "#117B79",
    "#117A78", "#117978", "#117878", "#117777", "#117677", "#127676", "#127576", "#127476", "#137375", "#137275", "#137174", "#147074",
    "#146F73", "#146E73", "#156D73", "#156C72", "#166B72", "#166A71", "#166971", "#176870", "#176770", "#17666F", "#18656F", "#18656E",
    "#18646E", "#19636D", "#19626D", "#19616C", "#19606C", "#1A5F6B", "#1A5E6B", "#1A5D6A", "#1A5C6A", "#1A5B69", "#1B5A68", "#1B5968",
    "#1B5867", "#1B5867", "#1B5766", "#1B5666", "#1C5565", "#1C5465", "#1C5364", "#1C5263", "#1C5163", "#1C5062", "#1C4F62", "#1C4E61",
    "#1C4D61", "#1C4C60", "#1C4C5F", "#1C4B5F", "#1C4A5E", "#1C495E", "#1C485D", "#1C475D", "#1C465C", "#1C455B", "#1C445B", "#1C435A",
    "#1C425A", "#1C4259", "#1C4158", "#1C4058", "#1B3F57", "#1B3E57", "#1B3D56", "#1B3C56", "#1B3B55", "#1B3A54", "#1B3954", "#1B3853",
    "#1A3753", "#1A3652", "#1A3651", "#1A3551", "#1A3450", "#1A3350", "#19324F", "#19314F", "#19304E", "#192F4D", "#192E4D", "#182D4C",
    "#182C4C", "#182B4B", "#182A4B", "#18294A", "#17284A", "#172749", "#172648", "#172548", "#172447", "#162347", "#162246", "#162146",
    "#162045", "#151F45", "#151E44", "#151D44"
]

function Base.show(io::IO, ::MIME"text/plain", rc::ResidueContributions)
    println(io, StyledStrings.styled"{bold:         Residue Contributions}")
    m = rc.residue_contributions
    min_val, max_val = extrema(m)
    if min_val == max_val
        max_val = min_val + 1.0
    end
    dstride = max(1, size(m, 1) ÷ 9 + 1)
    rstride = max(1, size(m, 2) ÷ 79 + 1)
    map = ""
    xlabel = false
    for d in size(m, 1):-dstride:1
        map *= if !xlabel && d < (size(m, 1) + dstride) ÷ 2
            xlabel = true
            " d  "
        else
            "    "
        end
        map *= "$(round(rc.d[d], digits=2)) "
        for res in 1:rstride:size(m, 2)
            cbin = tempo_color_scheme[
                    round(Int, 255 * (rc.residue_contributions[d, res] - min_val) / (max_val - min_val) + 1) 
                ]
            map *= StyledStrings.styled"{(fg=$cbin):█}"
        end
        map *= '\n'
    end
    map *= "         "
    for i in 1:rstride*8:length(rc.xticks[1])
        tick = "$(PDBTools.oneletter(rc.xticks[2][i][1:3]))$(rc.xticks[2][i][4:end])"
        tick *= repeat(" ", 8 - length(tick))
        map *= tick
    end
    map *= '\n'
    print(io, map)
end

function _custom_group_error_for_ResidueContributions()
    throw(ArgumentError("""\n

        ResidueContribution cannot be run if the MDDFs were computed with custom groups.

        This is because it requires the contribution of each independent atom to accumulate
        the residue contributions. The second argument of the function must be the vector of 
        atoms of the solute, not a group selection.

        Please read the documention, by typing: `? ResidueContribution`

    """))
end

# Custom error message for common mistaken call
ResidueContributions(result, g::Union{SoluteGroup,SolventGroup}, args...; kwargs...) =
    _custom_group_error_for_ResidueContributions()

@testitem "ResidueContribution" begin
    using ComplexMixtures
    using PDBTools
    using ComplexMixtures.Testing: data_dir, pdbfile

    atoms = readPDB(pdbfile)
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    water = AtomSelection(select(atoms, "water"), natomspermol=3)
    traj = Trajectory("$data_dir/NAMD/trajectory.dcd", protein, water)
    options = Options(;
        seed=1,
        stride=1,
        StableRNG=true,
        nthreads=2,
        silent=true
    )
    result = mddf(traj, options)

    # essential properties
    rc = ResidueContributions(result, select(atoms, "protein"))
    @test length(rc.d) == 101
    @test length.(rc.xticks) == (104, 104)
    @test size(rc.residue_contributions) == (101, 104)
    rc = ResidueContributions(result, select(atoms, "protein"); dmin=0.0, dmax=10.0)
    @test contributions(result, SoluteGroup(select(atoms, "protein and resnum 1"))) ≈ rc.residue_contributions[:, 1]
    @test contributions(result, SoluteGroup(select(atoms, "protein and resnum 104"))) ≈ rc.residue_contributions[:, 104]
    @test first(rc.d) == first(result.d)
    @test last(rc.d) == last(result.d)
    rcc = ResidueContributions(result, select(atoms, "protein"); dmin=0.0, dmax=10.0, type=:coordination_number)
    @test length(rcc.d) == 500
    @test length.(rcc.xticks) == (104, 104)
    @test size(rcc.residue_contributions) == (500, 104)
    @test contributions(result, SoluteGroup(select(atoms, "protein and resnum 1")); type=:coordination_number) ≈
          rcc.residue_contributions[:, 1]
    @test contributions(result, SoluteGroup(select(atoms, "protein and resnum 104")); type=:coordination_number) ≈
          rcc.residue_contributions[:, 104]

    # arithmetic operations
    rc = ResidueContributions(result, select(atoms, "protein"))
    rc2 = ResidueContributions(result, select(atoms, "protein"))
    @test all(rc.residue_contributions .- rc2.residue_contributions .< 1e-10)
    @test all(rc.residue_contributions + rc2.residue_contributions .≈ 2 .* rc.residue_contributions)
    rdiv = rc / rc2
    @test all(x -> isapprox(x, 1.0), filter(>(0.5), rdiv.residue_contributions))
    @test all(<(1.e-10), filter(<(0.5), rdiv.residue_contributions))
    rmul = rc * rc2
    @test rmul.residue_contributions ≈ rc.residue_contributions .^ 2

    # Error messages
    rc2 = ResidueContributions(result, select(atoms, "protein and residue < 50"))
    @test_throws ArgumentError rc - rc2
    rc2 = ResidueContributions(result, select(atoms, "protein"); dmax=3.0)
    @test_throws ArgumentError rc - rc2
    @test_throws ArgumentError ResidueContributions(result, select(atoms, "protein and resname XXX"))
    @test_throws ArgumentError ResidueContributions(result, SoluteGroup("ALA"))

    acidic_residues = select(atoms, "protein and acidic")
    basic_residues = select(atoms, "protein and basic")
    protein = AtomSelection(select(atoms, "protein"), nmols=1,
        group_atom_indices=[index.(acidic_residues), index.(basic_residues)],
        group_names=["acidic residues", "basic residues"]
    )
    traj = Trajectory("$data_dir/NAMD/trajectory.dcd", protein, water)
    result = mddf(traj, Options(lastframe=2))
    @test_throws ArgumentError ResidueContributions(result, select(atoms, "protein"))

end
