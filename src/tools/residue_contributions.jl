# These are type parameters necessary for indexing and iterating properly when a single
# residue contribution is selected, to iterate ver the contributions of such residue.
struct SingleResidueContribution end
struct MultipleResidueContribution end

"""
    ResidueContributions (data structure)

## Constructor function: 

    ResidueContributions(
        results::Result, atoms::AbstractVector{<:PDBTools.Atom};
        dmin=1.5, dmax=3.5,
        type=:mddf,
    )

Compute the residue contributions to the solute-solvent pair distribution function.
The function returns a `ResidueContributions` object that can be used to plot the residue contributions,
or to perform arithmetic operations with other `ResidueContributions` objects.

# Arguments

- `results::Result`: The result object obtained from the `mddf` function.
- `atoms::AbstractVector{<:PDBTools.Atom}`: The vector of atoms of the solute, or a part of it.

# Optional arguments

- `dmin::Float64`: The minimum distance to consider. Default is `1.5`.
- `dmax::Float64`: The maximum distance to consider. Default is `3.5`.
- `type::Symbol`: The type of the pair distribution function (`:mddf`, `:md_count`, or `:coordination_number`). Default is `:mddf`.
- `silent::Bool`: If `true`, the progress bar is not shown. Default is `false`.

A structure of type `ResultContributions` can be used to plot the residue contributions to the solute-solvent pair distribution function,
using the `Plots.contourf` function, and to perform arithmetic operations with other `ResidueContributions` objects, 
multiplying or dividing by a scalar, and slicing (see examples below). 

A ResidueContributions object can be saved to a JSON file using the `save` function, and loaded back using the `load` function.

# Examples

## Constructing a ResidueContributions object

```jldoctest
julia> using ComplexMixtures, PDBTools

julia> using ComplexMixtures.Testing: data_dir; ComplexMixtures._testing_show_method[] = true; # testing mode

julia> atoms = readPDB(data_dir*"/NAMD/Protein_in_Glycerol/system.pdb");

julia> results = load(data_dir*"/NAMD/Protein_in_Glycerol/protein_glyc.json");

julia> rc = ResidueContributions(results, select(atoms, "protein"); silent=true)

          Residue Contributions - 274 residues.
     3.51 █     █      █     █            █
     3.27 █              █   █
     3.03 █     █    █       █            █       █       █                █
     2.79 █    ██    █ █ █   █            █      ██          █        █    █
 d   2.55 █ █  ██    █ █ █   █            ██     ██ █  █  █  █  █     █    █
     2.31 █ █  ██    █ ███   █    ██      ██     ██ █  ██ █  ██ █     █    █
     2.07 █ █   █  █ █████   █    ██      ██     ██ █  ██ █  █ ██    █     █
     1.83 █   █ █  █ █████   █    ██      █      ██ █  ██    █ ██     █    █
     1.59
         A1      T33     T66     S98     S130    T162    A194    H226    G258     


julia> save("residue_contributions.json", rc)
"ResidueContributions saved in JSON file: residue_contributions.json"

julia> rc = load("residue_contributions.json", ResidueContributions);
```

## Plotting 

```julia
using ComplexMixtures, PDBTools, Plots
...
result = mddf(trajectory_file, solute, solvent, options)
rc = ResidueContributions(result, select(atoms, "protein"))
contourf(rc) # plots a contour map
```

## Slicing

Slicing, or indexing, the residue contributions returns a new `ResidueContributions` object with the selected residues:

```julia
using ComplexMixtures, PDBTools, Plots
...
result = mddf(trajectory_file, solute, solvent, options)
rc = ResidueContributions(result, select(atoms, "protein"))
rc_7 = rc[7] # contributions of residue 7
rc_range = rc[10:50] # slice the residue contributions
contourf(rc_range) # plots a contour map of the selected residues
```

## Single-residue contributions

When the contributions of a single residue are computed, or a single-residue contribution is retrieved from
a `ResidueContributions` object, the indexing and iteration over that object occurs over the contributions of that residue:

```julia
using ComplexMixtures, PDBTools, Plots
...
result = mddf(trajectory_file, solute, solvent, options)
rc = ResidueContributions(result, select(atoms, "protein"))
rc7 = rc[7] # contributions of residue 7
# iterate over the contributions of residue 7
rc7[1] # contribution of the first distance
rc7[end] # contribution of the last distance
```

This is particular useful to retrieve the contributions from all residues at a given distance:

```julia
rc = ResidueContributions(result, select(atoms, "protein"))
rc_last_distance = [ r[end] for r in rc ]
```

## Arithmetic operations

```julia
using ComplexMixtures, PDBTools, Plots
...
# first simulation (for example, low temperature):
result1 = mddf(trajectory_file1, solute, solvent, options)
rc1 = ResidueContributions(result1, select(atoms, "protein"))
# second simulation (for example, high temperature):
result2 = mddf(trajectory_file2, solute, solvent, options)
rc2 = ResidueContributions(result2, select(atoms, "protein"))
# difference of the residue contributions between the two simulations:
rc_diff = rc2 - rc1
contourf(rc_diff) # plots a contour map of the difference
# multiply and divide by scalars
rc3 = 3 * rc1
rc4 = rc2 / 2
```

!!! compat
    Slicing, indexing, and multiplication and divison by scalars were introduces in v2.7.0.
    Saving and loading was introduced in v2.8.0.

"""
@kwdef struct ResidueContributions{N<:Union{SingleResidueContribution, MultipleResidueContribution}}
    Version::VersionNumber = pkgversion(@__MODULE__)
    d::Vector{Float64}
    residue_contributions::Vector{Vector{Float64}}
    resnums::Vector{Int}
    xticks::Tuple{Vector{Int},Vector{String}}
end

function ResidueContributions(
    results::Result, atoms::AbstractVector{<:PDBTools.Atom};
    dmin=1.5, dmax=3.5,
    type=:mddf,
    silent=false,
    nthreads=Threads.nthreads(),
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
    rescontrib = [zeros(length(results.d)) for _ in 1:length(residues)]

    # Each column is then filled up with the contributions of each residue
    silent || (p = Progress(length(residues); dt=1))
    Threads.@threads for (ichunk, residue_inds) in enumerate(ChunkSplitters.index_chunks(residues; n=nthreads))
        _warn_zero_md_count = ichunk == 1 ? (!silent) : false
        for ires in residue_inds
            rescontrib[ires] .= contributions(results, SoluteGroup(residues[ires]); type, _warn_zero_md_count)
            _warn_zero_md_count = false
            silent || next!(p)
        end
    end

    # Plot only for distances within 1.5 and 3.5:
    # find the indexes of the distances that are within the range
    idmin = findfirst(d -> d > dmin, results.d)
    idmax = findfirst(d -> d > dmax, results.d)
    isnothing(idmin) && (idmin = 1)
    isnothing(idmax) && (idmax = length(results.d))

    # Obtain pretty labels for the residues in the x-axis (using PDBTools)
    resnums = PDBTools.resnum.(residues)
    xticks = PDBTools.residue_ticks(atoms;
        first=first(resnums),
        last=last(resnums),
        oneletter=false,
        serial=false,
    )

    N = length(resnums) == 1 ? SingleResidueContribution : MultipleResidueContribution
    return ResidueContributions{N}(;
        d=results.d[idmin:idmax],
        residue_contributions=[rc[idmin:idmax] for rc in rescontrib],
        resnums=resnums,
        xticks=xticks,
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

import Base: ==
==(rc1::ResidueContributions, rc2::ResidueContributions) =
    rc1.d == rc2.d && rc1.xticks == rc2.xticks && rc1.resnums == rc2.resnums &&
    rc1.residue_contributions == rc2.residue_contributions

function Base.copy(rc::RCType) where {RCType<:ResidueContributions}
    return RCType(
        d=copy(rc.d),
        residue_contributions=[copy(v) for v in rc.residue_contributions],
        resnums=copy(rc.resnums),
        xticks=(copy(rc.xticks[1]), copy(rc.xticks[2])),
    )
end


# 
# Indexing and iteration of MultipleResidueContribution objects: iterate over residues
#
function Base.getindex(rc::ResidueContributions{MultipleResidueContribution}, r::AbstractRange)
    N = length(r) == 1 ? SingleResidueContribution : MultipleResidueContribution
    return ResidueContributions{N}(
        d=rc.d,
        residue_contributions=[rc.residue_contributions[i] for i in r],
        resnums=rc.resnums[r],
        xticks=(rc.xticks[1][r], rc.xticks[2][r]),
    )
end
Base.length(rc::ResidueContributions{MultipleResidueContribution}) = length(rc.residue_contributions)
Base.getindex(rc::ResidueContributions{MultipleResidueContribution}, i) = rc[i:i]
Base.firstindex(rc::ResidueContributions{MultipleResidueContribution}) = 1
Base.lastindex(rc::ResidueContributions{MultipleResidueContribution}) = length(rc.residue_contributions)
function Base.iterate(rc::ResidueContributions{MultipleResidueContribution}, state=firstindex(rc))
    if state > length(rc)
        return nothing
    end
    return rc[state], state + 1
end

#
# Indexing and iteration of SingleResidueContribution objects: iterate over the contributions of a single residue
#
Base.length(rc::ResidueContributions{SingleResidueContribution}) = length(first(rc.residue_contributions))
Base.getindex(rc::ResidueContributions{SingleResidueContribution}, i) = first(rc.residue_contributions)[i]
Base.firstindex(rc::ResidueContributions{SingleResidueContribution}) = firstindex(first(rc.residue_contributions)) 
Base.lastindex(rc::ResidueContributions{SingleResidueContribution}) = length(first(rc.residue_contributions))
function Base.iterate(rc::ResidueContributions{SingleResidueContribution}, state=firstindex(first(rc.residue_contributions)))
    if state > length(first(rc.residue_contributions))
        return nothing
    end
    return first(rc.residue_contributions)[state], state + 1
end

import Base: -, +, /, *
function -(rc1::ResidueContributions, rc2::ResidueContributions)
    _check_identity_of_residues(rc1, rc2)
    rc_new = copy(rc1)
    for i in eachindex(rc_new.residue_contributions)
        @. rc_new.residue_contributions[i] = rc1.residue_contributions[i] - rc2.residue_contributions[i]
    end
    return rc_new
end
function +(rc1::ResidueContributions, rc2::ResidueContributions)
    _check_identity_of_residues(rc1, rc2)
    rc_new = copy(rc1)
    for i in eachindex(rc_new.residue_contributions)
        @. rc_new.residue_contributions[i] = rc1.residue_contributions[i] + rc2.residue_contributions[i]
    end
    return rc_new
end
function /(rc1::ResidueContributions, rc2::ResidueContributions)
    _check_identity_of_residues(rc1, rc2)
    if any(==(0), rc2.residue_contributions)
        @warn begin
            """\n
            Division by zero detected. These elements will be set to zero.

            """
        end _file = nothing _line = nothing
    end
    rc_new = copy(rc1)
    for i in eachindex(rc_new.residue_contributions)
        @. rc_new.residue_contributions[i] = ifelse(rc2.residue_contributions[i] == 0, 0.0, rc1.residue_contributions[i] / rc2.residue_contributions[i])
    end
    return rc_new
end
function *(rc1::ResidueContributions, rc2::ResidueContributions)
    _check_identity_of_residues(rc1, rc2)
    rc_new = copy(rc1)
    for i in eachindex(rc_new.residue_contributions)
        @. rc_new.residue_contributions[i] = rc1.residue_contributions[i] * rc2.residue_contributions[i]
    end
    return rc_new
end

#
# Arithmetic operations with scalars
#
function *(rc::ResidueContributions, x::Real)
    rc2 = copy(rc)
    for v in rc2.residue_contributions
        v .*= x
    end
    return rc2
end
*(x::Real, rc::ResidueContributions) = rc * x
/(rc::ResidueContributions, x::Real) = rc * inv(x)

const _colorscales = Dict{Symbol,Vector{Int}}(
    :tempo => [231, 194, 157, 120, 083, 046, 040, 034, 028, 022],
    :bwr => [017, 018, 019, 020, 021, 063, 105, 147, 189, 231, 224, 217, 210, 203, 196, 160, 124, 088, 052],
)

const _testing_show_method = Ref(false)

function _set_clims_and_colorscale!(rc::ResidueContributions; clims=nothing, colorscale=nothing)
    if isnothing(clims)
        minval = +Inf
        maxval = -Inf
        for v in rc.residue_contributions
            mi, ma = extrema(v)
            minval = min(minval, mi)
            maxval = max(maxval, ma)
        end
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
    if _testing_show_method[]
        print(io, "\n          Residue Contributions - $(length(rc)) residues.\n")
    else
        printstyled(io, "\n          Residue Contributions - $(length(rc)) residues.\n", bold=true)
    end
    m = rc.residue_contributions
    clims, colorscale = _set_clims_and_colorscale!(rc)
    colors = _colorscales[colorscale]
    ncolors = length(colors)
    dstride = max(1, length(rc.d) ÷ 9 + 1)
    rstride = max(1, length(m) ÷ 79 + 1)
    print(io, "")
    xlabel = false
    crange = clims[2] - clims[1]
    for d in length(rc.d):-dstride:1
        print(io, if !xlabel && d < (length(rc.d) + dstride) ÷ 2
            xlabel = true
            " d  "
        else
            "    "
        end)
        print(io, @sprintf("%5.2f ", rc.d[d]))
        for res in 1:rstride:length(m)
            cval = m[res][d]
            if _testing_show_method[]
                print(io, cval > 0.05 * clims[2] ? "█" : " ")
            else
                cbin = colors[round(Int, 1 + (ncolors - 1) * (cval - clims[1]) / crange)]
                printstyled(io, "█", color=cbin)
            end
        end
        println(io)
    end
    print(io, "         ")
    for i in 1:rstride*8:length(rc.xticks[1])
        tick = "$(PDBTools.oneletter(rc.xticks[2][i][1:3]))$(rc.resnums[i])"
        tick *= repeat(" ", 8 - length(tick))
        print(io, tick)
    end
    println(io)
end

function _custom_group_error_for_ResidueContributions()
    throw(ArgumentError("""\n

        The call to ResidueContributions with a SoluteGroup or SolventGroup objects as arguments is not valid. 

        You must provide *the atoms* of the solute as the second argument of the function,
        as a Vector{<:PDBTools.Atom} (a vector of Atom objects generated by PDBTools), for example:

        rc = ResidueContributions(result, PDBTools.select(atoms, "protein"))

        Please read the documention, by typing: ? ResidueContribution

        The computation of contributions of each residue requires independent atom contributions to 
        be stored, and thus custom groups are not supported. 

    """))
end

# Custom error message for common mistaken call
ResidueContributions(result, g::Union{SoluteGroup,SolventGroup}, args...; kwargs...) =
    _custom_group_error_for_ResidueContributions()

"""
    save(filename::String, rc::ResidueContributions)

Save the `ResidueContributions` object to a JSON file.

```julia
using ComplexMixtures
rc = ResidueContributions(resutls, SoluteGroup(protein))
save("residue_contributions.json", rc)
rc = load("residue_contributions.json", ResidueContributions)
```

"""
function save(filename::String, rc::ResidueContributions)
    filename = expanduser(filename)
    open(filename, "w") do f
        JSON3.write(f, rc)
    end
    return "ResidueContributions saved in JSON file: $filename"
end

"""
    load(filename::String, ResidueContributions)

Function to load the residue contributions saved into a JSON file into the `ResidueContributions` data structure.

## Example

```julia
using ComplexMixtures
rc = ResidueContributions(resutls, SoluteGroup(protein))
save("residue_contributions.json", rc)
rc = load("residue_contributions.json", ResidueContributions)
```

"""
function load(filename::String, ::Type{ResidueContributions})
    filename = expanduser(filename)
    _check_version(filename)
    rc = try
        open(filename, "r") do io
            JSON3.read(io, ResidueContributions{MultipleResidueContribution})
        end
    catch
        throw(ArgumentError("""\n
            The file $filename does not appear to contain a ResidueContributions object.

        """))
    end
    # Converto to SingleResidueContribution if a single contribution was read
    return if length(rc) == 1
        rc[1]
    else
        rc
    end
end

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
    @test rc.resnums == 1:104
    @test all(length(v) == 101 for v in rc.residue_contributions)
    @test length(rc.residue_contributions) == 104
    rc = ResidueContributions(result, select(atoms, "protein"); dmin=0.0, dmax=10.0)
    @test contributions(result, SoluteGroup(select(atoms, "protein and resnum 1"))) ≈ rc.residue_contributions[1]
    @test contributions(result, SoluteGroup(select(atoms, "protein and resnum 104"))) ≈ rc.residue_contributions[104]
    @test first(rc.d) == first(result.d)
    @test last(rc.d) == last(result.d)
    rcc = ResidueContributions(result, select(atoms, "protein"); dmin=0.0, dmax=10.0, type=:coordination_number)
    @test length(rcc.d) == 500
    @test length.(rcc.xticks) == (104, 104)
    @test all(length(v) == 500 for v in rcc.residue_contributions)
    @test length(rcc.residue_contributions) == 104
    @test contributions(result, SoluteGroup(select(atoms, "protein and resnum 1")); type=:coordination_number) ≈
          rcc.residue_contributions[1]
    @test contributions(result, SoluteGroup(select(atoms, "protein and resnum 104")); type=:coordination_number) ≈
          rcc.residue_contributions[104]

    # save and load
    tmpfile = tempname() * ".json"
    save(tmpfile, rc)
    rc_load = load(tmpfile, ResidueContributions)
    @test rc == rc_load
    tmpfile = tempname()
    open(tmpfile, "w") do io
        println(io, """{"Version":"$(pkgversion(@__MODULE__))"}""")
    end
    @test_throws ArgumentError load(tmpfile, ResidueContributions)

    # version issues
    N = ComplexMixtures.MultipleResidueContribution
    rc_future = ResidueContributions{N}(Version=v"1000.0.0", d=rc.d, residue_contributions=rc.residue_contributions, resnums=rc.resnums, xticks=rc.xticks)
    save(tmpfile, rc_future)
    @test_throws ArgumentError load(tmpfile, ResidueContributions)
    rc_past = ResidueContributions{N}(Version=v"1.0.0", d=rc.d, residue_contributions=rc.residue_contributions, resnums=rc.resnums, xticks=rc.xticks)
    save(tmpfile, rc_past)
    @test_throws ArgumentError load(tmpfile, ResidueContributions)

    # indexing
    rc = ResidueContributions(result, select(atoms, "protein"))
    rc1 = rc[1]
    @test rc1.resnums == [1]
    @test rc1 == ResidueContributions(result, select(atoms, "protein and resnum 1"))
    rc2 = rc[2:10]
    @test rc2.resnums == 2:10
    @test rc2 == ResidueContributions(result, select(atoms, "protein and resnum > 1 and resnum < 11"))
    @test firstindex(rc) == 1
    @test lastindex(rc) == 104
    @test length(rc) == 104
    @test firstindex(rc1) == 1
    @test lastindex(rc1) == 101
    @test length(rc1) == 101

    # Test iterators
    @test count(true for r in rc) == 104
    @test count(true for r in rc1) == 101
    @test [ r[end] for r in rc ] == [ r[end] for r in rc.residue_contributions ]

    # empty plot (just test if the show function does not throw an error)
    rc2 = copy(rc)
    for v in rc2.residue_contributions
        v .= 0.0
    end
    @test show(IOBuffer(), MIME"text/plain"(), rc2) === nothing

    # arithmetic operations
    rc = ResidueContributions(result, select(atoms, "protein"))
    rcminus = rc - rc
    @test all(sum(rcminus.residue_contributions) .< 1e-10)
    rcplus = rc + rc
    @test all((rcplus.residue_contributions[i] ≈ 2 * rc.residue_contributions[i]) for i in 1:length(rc.residue_contributions))
    rdiv = rc / rc
    @test all(all(x -> isapprox(x, 1.0), filter(>(0.0), v)) for v in rdiv.residue_contributions)
    @test all(all(<(1e-10), filter(<(0.0), v)) for v in rdiv.residue_contributions)
    rmul = rc * rc
    @test all(rmul.residue_contributions[i] ≈ rc.residue_contributions[i] .^ 2 for i in 1:length(rc.residue_contributions))
    rc2 = 2 * rc
    @test all(rc2.residue_contributions[i] == 2 .* rc.residue_contributions[i] for i in 1:length(rc.residue_contributions))
    rc2 = rc * 2
    @test all(rc2.residue_contributions[i] == 2 .* rc.residue_contributions[i] for i in 1:length(rc.residue_contributions))
    rc2 = rc / 2
    @test all(rc2.residue_contributions[i] == rc.residue_contributions[i] ./ 2 for i in 1:length(rc.residue_contributions))

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

    # Solute with more than one molecule
    glicines = select(atoms, "protein and resname GLY")
    solute = AtomSelection(glicines; natomspermol=7)
    water = AtomSelection(select(atoms, "water"), natomspermol=3)
    traj = Trajectory("$data_dir/NAMD/trajectory.dcd", solute, water)
    options = Options(;
        seed=1,
        stride=1,
        StableRNG=true,
        nthreads=2,
        silent=true
    )
    result = mddf(traj, options)
    rc = ResidueContributions(result, glicines)
    # This might actually be changed in the future, depending on what one wants. Maybe just throw an error.
    @test all(rc.residue_contributions[i] ≈ rc.residue_contributions[1] for i in 1:length(rc.residue_contributions))

end
