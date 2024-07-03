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
Base.copy(rc::ResidueContributions) = 
    ResidueContributions(copy(rc.d), (copy(rc.xticks[1]), copy(rc.xticks[2])), copy(rc.residue_contributions))

const _colorscales = Dict{Symbol,Vector{Int}}(
    :tempo => [231, 194, 157, 120, 083, 046, 040, 034, 028, 022 ],
    :bwr => [017, 018, 019, 020, 021, 063, 105, 147, 189, 231, 224, 217, 210, 203, 196, 160, 124, 088, 052 ],
)

function _set_clims_and_colorscale!(rc::ResidueContributions; clims=nothing, colorscale=nothing)
    if isnothing(clims)
        minval, maxval = extrema(rc.residue_contributions)
        if minval == maxval
            maxval = minval + 1.0
        end
        if minval >= 0
            clims = (0, maxval)
            isnothing(colorscale) && (colorscale = :tempo)
        elseif maxval <= 0
            clims = (minval, 0)
            isnothing(colorscale) && (colorscale = :tempo)
        else
            mval = max(abs(minval), maxval)
            clims = (-mval, mval)
            isnothing(colorscale) && (colorscale = :bwr)
        end
    else
        colorscale = isnothing(colorscale) ? :bwr : colorscale
    end
    return clims, colorscale
end

function Base.show(io::IO, ::MIME"text/plain", rc::ResidueContributions)
    printstyled(io, "\n          Residue Contributions\n", bold=true)
    m = rc.residue_contributions
    clims, colorscale = _set_clims_and_colorscale!(rc)
    colors = _colorscales[colorscale]
    ncolors = length(colors)
    dstride = max(1, size(m, 1) ÷ 9 + 1)
    rstride = max(1, size(m, 2) ÷ 79 + 1)
    print(io,"")
    xlabel = false
    crange = clims[2] - clims[1]
    for d in size(m, 1):-dstride:1
        print(io, if !xlabel && d < (size(m, 1) + dstride) ÷ 2
            xlabel = true
            " d  "
        else
            "    "
        end)
        print(io, @sprintf("%5.2f ", rc.d[d]))
        for res in 1:rstride:size(m, 2)
            cval = rc.residue_contributions[d, res]
            cbin = colors[round(Int, 1 + (ncolors - 1) * (cval - clims[1]) / crange)]
            printstyled(io, "█", color=cbin)
        end
        println(io)
    end
    print(io,"         ")
    for i in 1:rstride*8:length(rc.xticks[1])
        tick = "$(PDBTools.oneletter(rc.xticks[2][i][1:3]))$(rc.xticks[2][i][4:end])"
        tick *= repeat(" ", 8 - length(tick))
        print(io,tick)
    end
    println(io)    
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

    # copy structure
    rc2 = copy(rc)
    @test rc2.d == rc.d
    @test rc2.residue_contributions == rc.residue_contributions
    @test rc2.xticks[1] == rc.xticks[1]
    @test rc2.xticks[2] == rc.xticks[2]

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
