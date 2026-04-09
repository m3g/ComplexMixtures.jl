"""

$(TYPEDEF)

Structure to contain the density values obtained from the calculation.

$(TYPEDFIELDS)

"""
@kwdef mutable struct Density
    solute::Float64 = 0.0
    solvent::Float64 = 0.0
    solvent_bulk::Float64 = 0.0
end

function Base.show(io::IO, d::Density)
    println(io, "Density of solute: ", d.solute)
    println(io, "Density of solvent: ", d.solvent)
    print(io, "Density of solvent in bulk: ", d.solvent_bulk)
end

#=
    setbin(d,step)

Function that sets to which histogram bin a data point pertains simple, but important to keep consistency over all calls.

=#
setbin(d, step) = max(1, ceil(Int, d / step))

#=

Structures to contain the volumes obtained from calculations.

=#
@kwdef mutable struct Volume
    total::Float64
    bulk::Float64
    domain::Float64
    shell::Vector{Float64}
end
Volume(nbins::Int) = Volume(0.0, 0.0, 0.0, zeros(Float64, nbins))

function Base.show(io::IO, v::Volume)
    println(io, "Total volume: ", v.total)
    println(io, "Bulk volume: ", v.bulk)
    println(io, "Domain volume: ", v.domain)
    print(io, "Shell volumes: ", v.shell)
end

"""

$(TYPEDEF)

Structure to contain the results of the MDDF calculation.

$(TYPEDFIELDS)

The Result{Vector{Float64}} parametric type is necessary only for reading the JSON3 saved file. 

"""
@kwdef mutable struct Result

    # ComplexMixtures version that generated this results
    Version::VersionNumber = pkgversion(@__MODULE__)

    # Histogram properties
    nbins::Int
    dbulk::Float64
    cutoff::Float64
    d::Vector{Float64} = zeros(nbins)

    # Data to compute the MDDF distribution and corresponding KB integral
    md_count::Vector{Float64} = zeros(nbins)
    md_count_random::Vector{Float64} = zeros(nbins)
    coordination_number::Vector{Float64} = zeros(nbins)
    coordination_number_random::Vector{Float64} = zeros(nbins)

    # The resulting MDDF and kb integrals
    mddf::Vector{Float64} = zeros(nbins)
    kb::Vector{Float64} = zeros(nbins)

    # Properties of the solute and solvent selections
    autocorrelation::Bool
    solute::AtomSelection
    solvent::AtomSelection

    # Group (atomic type by default) contributions to 
    # the coordination number counts. These are used to
    # compute group contributions to the MDDFs and KBIs
    # Note: These could be Matrix{Float64}, but for the convenience
    # of using JSON3, we use Vector{Vector{Float64}}, which is 
    # read directly.
    solute_group_count::Vector{Vector{Float64}}
    solvent_group_count::Vector{Vector{Float64}}
    solute_group_count_random::Vector{Vector{Float64}}
    solvent_group_count_random::Vector{Vector{Float64}}

    # Data to compute a RDF and the KB integral from this count
    rdf_count::Vector{Float64} = zeros(nbins)
    rdf_count_random::Vector{Float64} = zeros(nbins)
    sum_rdf_count::Vector{Float64} = zeros(nbins)
    sum_rdf_count_random::Vector{Float64} = zeros(nbins)
    rdf::Vector{Float64} = zeros(nbins)
    kb_rdf::Vector{Float64} = zeros(nbins)

    # Overall densities and volumes
    density::Density = Density()
    volume::Volume = Volume(nbins)

    #
    # If multiple files, the following vector defines the weight of each file
    # in the final merged results. It is proportional to the number of frames
    # and weights provided for the frames in each file. 
    #
    files::Vector{TrajectoryFileOptions}
    weights::Vector{Float64}

end

#
# Initialize the data structure that is returned from the computation, and checks some
# input parameters for consistency
#
function Result(
    trajectory::Trajectory,
    options::Options;
    trajectory_data::TrajectoryMetaData=TrajectoryMetaData(trajectory, options),
    frame_weights::AbstractArray{<:Real}=Float64[]
)

    # Number of bins of the histogram
    nbins = setbin(options.cutoff, options.binstep)

    # If frame weights are provided, the length of the weights vector has to at least
    # of the of number of the last frame to be read, and the sum of the weights of the
    # frames to be considered must be greater than zero
    if !isempty(frame_weights)
        if count(>(1), size(frame_weights)) > 1
            throw(ArgumentError(chomp("""\n
                The frame_weights provided must be one-dimensional, but a $(size(frame_weights))-dimensional array was given.

            """)))
        end
        if length(frame_weights) < trajectory_data.lastframe_read
            throw(ArgumentError(chomp("""\n
            The length of the frame_weights vector provided must at least the number of frames to be read.

                Input given: length(frame_weights) = $(length(frame_weights))
                             last frame to be read: $(trajectory_data.lastframe_read)

            """)))
        end
        range_considered = options.firstframe:options.stride:trajectory_data.lastframe_read
        if sum(frame_weights[i] for i in range_considered) <= 0.0
            throw(ArgumentError("""\n
                The sum of the frame weights must be greater than zero for the frames that will be considered: $(range_considered)

            """))
        end
    else
        frame_weights = ones(trajectory_data.lastframe_read)
    end

    # Initialize group count arrays
    solute_group_count = [zeros(nbins) for _ in 1:trajectory_data.n_groups_solute]
    solvent_group_count = [zeros(nbins) for _ in 1:trajectory_data.n_groups_solvent]
    solute_group_count_random = [zeros(nbins) for _ in 1:trajectory_data.n_groups_solute]
    solvent_group_count_random = [zeros(nbins) for _ in 1:trajectory_data.n_groups_solvent]

    return Result(
        nbins=nbins,
        dbulk=options.dbulk,
        cutoff=options.cutoff,
        autocorrelation=isautocorrelation(trajectory),
        solute=trajectory.solute,
        solvent=trajectory.solvent,
        solute_group_count=solute_group_count,
        solvent_group_count=solvent_group_count,
        solute_group_count_random=solute_group_count_random,
        solvent_group_count_random=solvent_group_count_random,
        files=[
            TrajectoryFileOptions(
                filename=trajectory.filename,
                options=options,
                irefatom=trajectory_data.irefatom,
                lastframe_read=trajectory_data.lastframe_read,
                nframes_read=trajectory_data.nframes_read,
                frame_weights=[Float64(frame_weights[i]) for i in eachindex(frame_weights)],
            )
        ],
        weights=[1.0],
    )
end

@testitem "input: argument errors" begin
    using ComplexMixtures: mddf, Trajectory, Options, AtomSelection
    using ComplexMixtures: data_dir
    using PDBTools: read_pdb, select
    atoms = read_pdb("$data_dir/NAMD/structure.pdb")
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    water = AtomSelection(select(atoms, "water"), natomspermol=3)
    trajectory = Trajectory("$data_dir/NAMD/trajectory.dcd", tmao, water)
    @test_throws ArgumentError mddf(trajectory, Options(lastframe=100))
    @test_throws ArgumentError mddf(trajectory, Options(irefatom=1000))
    @test_throws ArgumentError mddf(trajectory, Options(lastframe=5), frame_weights=[1.0, 1.0, 1.0])
    options = Options(firstframe=2, lastframe=5, stride=2)
    @test_throws ArgumentError mddf(trajectory, options, frame_weights=[1.0, 0.0, 1.0, 0.0, 1.0])
end

#
# What to show at the REPL
#
Base.show(io::IO, R::Result) = show(io, overview(R))

#=
Type of calculation: auto or cross-correlations - have different sampling schemes

Choose self and cross versions

In both the self and (cross) non-self cases, the number of random samples is 
n_random_samples*solvent.nmols. However, in the non-self distribution
the sampling of solvent distances is proportional to the number of solute
molecules, thus the md count has to be averaged by solute.nmols. In the
case of the self-distribution, we compute n(n-1)/2 distances, and we will
divide this by n random samples, which is the sampling of the random
distribution. Therefore, we must weight the self-distance count by dividing
it by (n-1)/2, so that we have a count proportional to n as well, leading
to the correct weight relative to the random sample. 
=#
function set_samples(R::Result)
    if R.autocorrelation
        samples = (solvent_nmols=R.solvent.nmols - 1, random=R.files[1].options.n_random_samples)
    else
        samples = (solvent_nmols=R.solvent.nmols, random=R.files[1].options.n_random_samples)
    end
    return samples
end

# autocorrelation can be obtained from the comparison of solute and solvent indices
isautocorrelation(solute_indices::AbstractVector{<:Integer}, solvent_indices::AbstractVector{<:Integer}) =
    solute_indices == solvent_indices ? true : false
isautocorrelation(trajectory::Trajectory) =
    isautocorrelation(trajectory.solute.indices, trajectory.solvent.indices)

#
# Functions to compute volumes of shells
#
#=
    sphericalshellvolume(i,step)

Computes the volume of the spherical shell defined within [(i-1)*step,i*step].

=#
function sphericalshellvolume(i, step)
    rmin = (i - 1) * step
    return (4 * pi / 3) * ((rmin + step)^3 - rmin^3)
end

@testitem "sphericalshellvolume" begin
    using ComplexMixtures: sphericalshellvolume
    @test sphericalshellvolume(1, 1.0) ≈ 4 * pi / 3
    @test sphericalshellvolume(2, 1.0) ≈ 4 * pi / 3 * (8 - 1)
    @test sphericalshellvolume(3, 1.0) ≈ 4 * pi / 3 * (27 - 8)
end

#=
    shellradius(i,step)

Compute the point in which the radius comprises half of the volume of the shell.

=#
function shellradius(i, step)
    rmin = (i - 1) * step
    return (0.5 * ((rmin + step)^3 + rmin^3))^(1 / 3)
end

@testitem "shellradius" begin
    using ComplexMixtures: shellradius
    @test shellradius(1, 0.1) ≈ 0.07937005259840998
    @test shellradius(5, 0.3) ≈ 1.3664650373440481
end

#
# Function to compute the normalization of the weights, given the optional
# frame weights and the number of frames read, which is used in the case
# that frame weights were not provided
#
function sum_frame_weights(R::Result)
    Q = 0.0
    for file in R.files
        i = file.options.firstframe
        s = file.options.stride
        l = file.lastframe_read
        Q += sum(file.frame_weights[i] for i in i:s:l)
    end
    return Q
end

#=
    finalresults!(R::Result, options::Options; coordination_number_only)

Function that computes the final results of all the data computed by averaging according to the sampling of each type of data, and converts to common units.

Computes also the final distribution functions and KB integrals.

This function modified the values contained in the R data structure.

If `coordination_number_only` is true, do not print a warning if there are zero samples in the ideal-gas histogram bins.

=#
function finalresults!(R::Result, options::Options; coordination_number_only)
    R = if !coordination_number_only
        _mddf_final_results!(R, options)
    else
        _coordination_number_final_results!(R, options)
    end
    return R
end

function _mddf_final_results!(R::Result, options::Options)

    # Setup the distance vector
    R.d .= shellradius.(1:R.nbins, options.binstep)

    # Sampling scheme depending on the type of calculation
    samples = set_samples(R)

    # Normalization of number of frames: sum of weights for all frames read
    Q = sum_frame_weights(R)

    # Scale counters by number of samples and frames
    @. R.md_count = R.md_count / (R.solute.nmols * Q)
    @. R.solute_group_count = R.solute_group_count / (R.solute.nmols * Q)
    @. R.solute_group_count_random = R.solute_group_count_random / (samples.random * Q)
    if R.autocorrelation
        R.solvent_group_count .= R.solute_group_count
        R.solvent_group_count_random .= R.solute_group_count_random 
    else
        @. R.solvent_group_count = R.solvent_group_count / (samples.random * Q)
        @. R.solvent_group_count_random = R.solvent_group_count_random / (samples.random * Q)
    end
    @. R.md_count_random = R.md_count_random / (samples.random * Q)
    @. R.rdf_count = R.rdf_count / (R.solute.nmols * Q)
    @. R.rdf_count_random = R.rdf_count_random / (samples.random * Q)

    # Volume of each bin shell and of the solute domain
    R.volume.total = R.volume.total / Q
    @. R.volume.shell = R.volume.total * (R.rdf_count_random / samples.solvent_nmols)

    # Solute domain volume
    ibulk = setbin(R.dbulk + 0.5 * R.files[1].options.binstep, R.files[1].options.binstep)
    R.volume.domain = sum(@view(R.volume.shell[1:ibulk-1]))

    # Bulk volume and density properties: either the bulk is considered everything
    # that is not the domain, or the bulk is the region between d_bulk and cutoff,
    # if R.files[1].options.usecutoff is true (meaning that there is a cutoff different from
    # that of the bulk distance)
    if !R.files[1].options.usecutoff
        R.volume.bulk = R.volume.total - R.volume.domain
        n_solvent_in_bulk = samples.solvent_nmols - sum(R.rdf_count)
    else
        n_solvent_in_bulk = sum(@view(R.rdf_count[ibulk:R.nbins]))
        R.volume.bulk = sum(@view(R.volume.shell[ibulk:R.nbins]))
    end
    R.density.solvent = R.solvent.nmols / R.volume.total
    R.density.solute = R.solute.nmols / R.volume.total
    R.density.solvent_bulk = n_solvent_in_bulk / R.volume.bulk

    # Now that we know the the volume of the domain and the density of the solvent in the 
    # bulk region, we can rescale the random counts to take into account that the the ideal
    # gas distribution must have more molecules, with same bulk density, that the true
    # distribution, because we have to take into consieration the available volume which is
    # occupied by the solute
    density_fix = R.density.solvent_bulk / R.density.solvent
    return renormalize!(R, density_fix; silent=options.silent)
end

function renormalize!(R::Result, density_fix::Number; silent)
    R.md_count_random .*= density_fix 
    R.rdf_count_random .*= density_fix
    R.solute_group_count_random .*= density_fix
    R.solvent_group_count_random .*= density_fix

    # Computing the distribution functions and KB integrals, from the MDDF and from the RDF
    #
    warned_already = false
    cumsum!(R.coordination_number, R.md_count)
    cumsum!(R.coordination_number_random, R.md_count_random)
    for ibin = 1:R.nbins
        if R.md_count_random[ibin] > 0.0
            R.mddf[ibin] = R.md_count[ibin] / R.md_count_random[ibin]
        else
            if !warned_already && !silent
                @warn begin
                    """\n
                        Ideal-gas histogram bins with zero samples. 
                        Increase n_random_samples, number of trajectory frames, and/or bin size.

                    """
                end _file = nothing _line = nothing
                warned_already = true
            end
        end
        R.kb[ibin] =
            units.Angs3tocm3permol *
            (1 / R.density.solvent_bulk) *
            (R.coordination_number[ibin] - R.coordination_number_random[ibin])

        # For the RDF
        if R.rdf_count_random[ibin] > 0.0
            R.rdf[ibin] = R.rdf_count[ibin] / R.rdf_count_random[ibin]
        end
        if ibin == 1
            R.sum_rdf_count[ibin] = R.rdf_count[ibin]
            R.sum_rdf_count_random[ibin] = R.rdf_count_random[ibin]
        else
            R.sum_rdf_count[ibin] = R.sum_rdf_count[ibin-1] + R.rdf_count[ibin]
            R.sum_rdf_count_random[ibin] =
                R.sum_rdf_count_random[ibin-1] + R.rdf_count_random[ibin]
        end
        R.kb_rdf[ibin] =
            units.Angs3tocm3permol *
            (1 / R.density.solvent_bulk) *
            (R.sum_rdf_count[ibin] - R.sum_rdf_count_random[ibin])

    end
    return R
end

function _coordination_number_final_results!(R::Result, options::Options)

    if !options.silent
        @warn begin
            """\n
                coordination_number_only was set to true, so the MDDF and KB integrals were not computed. 
                (to remove this warning use `Options(silent=true)`)

            """
        end _file = nothing _line = nothing
    end

    # Setup the distance vector
    R.d .= shellradius.(1:R.nbins, options.binstep)

    # Normalization of number of frames: sum of weights for all frames read
    Q = sum_frame_weights(R)

    # Scale counters by number of samples and frames
    @. R.md_count = R.md_count / (R.solute.nmols * Q)
    @. R.solute_group_count = R.solute_group_count / (R.solute.nmols * Q)
    if R.autocorrelation
        R.solvent_group_count .= R.solute_group_count
    else
        @. R.solvent_group_count = R.solvent_group_count / (R.solute.nmols * Q)
    end
    @. R.rdf_count = R.rdf_count / (R.solute.nmols * Q)

    # Volume of each bin shell and of the solute domain
    R.volume.total = R.volume.total / Q

    R.density.solvent = R.solvent.nmols / R.volume.total
    R.density.solute = R.solute.nmols / R.volume.total

    # Coordination numbers
    cumsum!(R.coordination_number, R.md_count)
    cumsum!(R.sum_rdf_count, R.rdf_count)

    return R
end

@testitem "Result - empty" begin
    #using ComplexMixtures: Result, Trajectory, Options, AtomSelection
    using ComplexMixtures: data_dir, pdb_file_example
    using PDBTools: read_pdb, select, name
    atoms = read_pdb(pdb_file_example)
    protein = select(atoms, "protein")
    tmao = select(atoms, "resname TMAO")
    solute = AtomSelection(protein, nmols=1)
    solvent = AtomSelection(tmao, natomspermol=14)
    traj = Trajectory("$data_dir/NAMD/trajectory.dcd", solute, solvent)
    options = Options()
    # At this point we can only test an empty Result struct
    R = Result(traj, options)
    @test R.autocorrelation == false
    @test R.cutoff == 10.0
    @test length(R.d) == 500
    @test R.dbulk == 10.0
    @test (R.density.solute, R.density.solvent, R.density.solvent_bulk) == (0.0, 0.0, 0.0)
    @test normpath(R.files[1].filename) == normpath("$data_dir/NAMD/trajectory.dcd")
    @test R.files[1].irefatom == 1
    @test length(R.kb) == 500
    @test length(R.kb_rdf) == 500
    @test R.files[1].lastframe_read == 20
    @test length(R.md_count) == 500
    @test length(R.md_count_random) == 500
    @test length(R.mddf) == 500
    @test R.nbins == 500
    @test R.files[1].nframes_read == 20
    @test R.files[1].options == Options()
    @test length(R.rdf_count) == 500
    @test length(R.rdf_count_random) == 500
    @test R.solute == AtomSelection(collect(1:1463), nmols=1, natomspermol=1463, group_names=name.(protein))
    @test R.solvent == AtomSelection(collect(1479:4012), nmols=181, natomspermol=14, group_names=name.(tmao[1:14]))
    @test length(R.solute_group_count) == 1463
    @test all(length.(R.solute_group_count) .== 500)
    @test length(R.solvent_group_count) == 14
    @test all(length.(R.solvent_group_count) .== 500)
    @test length(R.coordination_number) == 500
    @test length(R.coordination_number_random) == 500
    @test length(R.sum_rdf_count_random) == 500
    @test (R.volume.total, R.volume.bulk, R.volume.domain, length(R.volume.shell)) == (0.0, 0.0, 0.0, 500)
    @test R.weights[1] == 1.0
end

#
# Functions to save the results to a file
#
# These definitions are probably not necessary for StructTypes >= 1.10.0
# 
StructTypes.StructType(::Type{AtomSelection}) = StructTypes.Struct()
StructTypes.StructType(::Type{Result}) = StructTypes.Struct()
StructTypes.StructType(::Type{Density}) = StructTypes.Struct()
StructTypes.StructType(::Type{Volume}) = StructTypes.Struct()
StructTypes.StructType(::Type{Options}) = StructTypes.Struct()
StructTypes.StructType(::Type{TrajectoryFileOptions}) = StructTypes.Struct()

"""
    save(filename::AbstractString, R::Result)

Function to write the result data structure to a json file.

"""
function save(filename::AbstractString, R::Result)
    filename = expanduser(filename)
    open(filename, "w") do f
        JSON3.write(f, R)
    end
    return "Results saved in JSON file: $filename"
end
# legacy order
save(R::Result, filename::AbstractString) = save(filename, R)

#
# This function tries to read a version number from a result.json
# file. This is used to adjust the reading to legacy formats of the
# output file.
#
function _get_version(filename)
    str = readuntil(expanduser(filename), ',')
    v = match(r"\"Version\":\"([^\"]*)\"", str)
    return isnothing(v) ? v"1.0.0" : VersionNumber(v[1])
end

"""
    load(filename::AbstractString, [::Type{Result}=Result])

Function to load the json saved results file into, by default, the `Result` data structure.
The second parameter is optional for loading `Result` objects.

## Example

```julia
using ComplexMixtures: load
R = load("results.json")
#or
R = load("results.json", Result)
```

"""
function load(filename::AbstractString, ::Type{Result})
    filename = expanduser(filename)
    _check_version(filename)
    R = try
        open(filename, "r") do io
            JSON3.read(io, Result)
        end
    catch
        throw(ArgumentError("""\n 
            The file $filename does not appear to contain a Result object.

        """))
    end
    return R
end
load(filename::AbstractString) = load(filename, Result)

@testitem "Result - load/save" begin
    using ComplexMixtures: load
    using ComplexMixtures: data_dir
    r1 = load("$data_dir/NAMD/protein_tmao.json")
    tmp = tempname()
    save(r1, tmp)
    r2 = load(tmp)
    @test r1 == r2
    r2 = load(tmp, Result)
    @test r1 == r2
    save(tmp, r1)
    r2 = load(tmp)
    @test r1 == r2
    # Test throwing an error incompatible versions of ComplexMixtures
    @test_throws ArgumentError load("$data_dir/wrong_version_jsons/too_new.json")
    @test_throws ArgumentError load("$data_dir/wrong_version_jsons/too_old.json")
    rm(tmp)
    tmpfile = tempname()
    open(tmpfile, "w") do io
        println(io, """{"Version":"$(pkgversion(@__MODULE__))"}""")
    end
    @test_throws ArgumentError load(tmpfile)
    rm(tmpfile)

    # Test load and save with substrings
    r1 = load(@view("$data_dir/NAMD/protein_tmao.json---"[1:end-3]))
    tmp = @view(tempname()[1:end-1])
    save(r1, tmp)
    r2 = load(tmp)
    @test r1 == r2
end

#=
    sum!(R1::Result, R2::Result)

Sum the counts of two Results structures, adding the result to the first structure as in R1 = R1 + R2.

=#
function sum!(R1::Result, R2::Result)

    @. R1.md_count += R2.md_count
    @. R1.md_count_random += R2.md_count_random

    for i in eachindex(R1.solute_group_count)
        @. R1.solute_group_count[i] += R2.solute_group_count[i]
        @. R1.solute_group_count_random[i] += R2.solute_group_count_random[i]
    end
    for i in eachindex(R1.solvent_group_count)
        @. R1.solvent_group_count[i] += R2.solvent_group_count[i]
        @. R1.solvent_group_count_random[i] += R2.solvent_group_count_random[i]
    end

    @. R1.rdf_count += R2.rdf_count
    @. R1.rdf_count_random += R2.rdf_count_random

    R1.volume.total += R2.volume.total

    return R1
end
#=
    title(R::Result, solute::AtomSelection, solvent::AtomSelection)
    title(R::Result, solute::AtomSelection, solvent::AtomSelection, nspawn::Int)

Print some information about the run.

=#
function title(R::Result, solute::AtomSelection, solvent::AtomSelection)
    print(
        """
        $(bars)
        Starting MDDF calculation:
        $(R.files[1].nframes_read) frames will be considered.
        Solute: $(atoms_str(natoms(solute))) belonging to $(mol_str(solute.nmols)).
        Solvent: $(atoms_str(natoms(natoms))) belonging to $(mol_str(solvent.nmols))
        """)
end
function title(R::Result, solute::AtomSelection, solvent::AtomSelection, nspawn::Int)
    print(
        """ 
        $(bars)
        Starting MDDF calculation in parallel:
        $(R.files[1].nframes_read) frames will be considered.
        Number of calculation threads: $(nspawn)
        Solute: $(atoms_str(natoms(solute))) belonging to $(mol_str(solute.nmols)).
        Solvent: $(atoms_str(natoms(solvent))) belonging to $(mol_str(solvent.nmols)).
        """)
end

#
# Print overview of the results in the REPL
#
#=

Structure that is used to dispatch the show of a overview.

=#
@kwdef mutable struct Overview
    R::Result
    domain_molar_volume::Float64 = 0.0
    density::Density = Density()
    solvent_molar_volume::Float64 = 0.0
    solvent_molar_volume_bulk::Float64 = 0.0
    solute_molar_volume::Float64 = 0.0
end

function _bulk_range_from_R(R)
    if R.dbulk == R.cutoff
        return ">= $(R.dbulk) Å"
    else
        return "$(R.dbulk) - $(R.cutoff) Å"
    end
end

function Base.show(io::IO, ov::Overview)
    println(
        io,
        """
 $bars
 MDDF Overview - ComplexMixtures - Version $(ov.R.Version)
 $bars

 Solvent properties:
 -------------------

 Simulation concentration: $(ov.density.solvent) mol L⁻¹
 Molar volume: $(ov.solvent_molar_volume) cm³ mol⁻¹

 Concentration in bulk: $(ov.density.solvent_bulk) mol L⁻¹
 Molar volume in bulk: $(ov.solvent_molar_volume_bulk) cm³ mol⁻¹

 Solute properties:
 ------------------

 Simulation Concentration: $(ov.density.solute) mol L⁻¹
 Estimated solute partial molar volume: $(ov.solute_molar_volume) cm³ mol⁻¹

 Bulk range: $(_bulk_range_from_R(ov.R))
 Molar volume of the solute domain: $(ov.domain_molar_volume) cm³ mol⁻¹

 Auto-correlation: $(ov.R.autocorrelation)

 Trajectory files and weights:
 """
    )
    for i in eachindex(ov.R.files)
        println(io, "   $(normpath(ov.R.files[i].filename)) - w = $(ov.R.weights[i])")
    end
    ifar = trunc(Int, ov.R.nbins - 1.0 / ov.R.files[1].options.binstep)
    long_range_mean = mean(ov.R.mddf[ifar:ov.R.nbins])
    long_range_std = std(ov.R.mddf[ifar:ov.R.nbins])
    long_range_mean_rdf = mean(ov.R.rdf[ifar:ov.R.nbins])
    long_range_std_rdf = std(ov.R.rdf[ifar:ov.R.nbins])
    print(
        io,
        """

   Long range MDDF mean (expected 1.0): $long_range_mean ± $long_range_std
   Long range RDF mean (expected 1.0): $long_range_mean_rdf ± $long_range_std_rdf

   $bars"""
    )
end

"""
    overview(R::Result)

Function that outputs the volumes and densities in the most natural units.
"""
function overview(R::Result)

    ov = Overview(R=R)

    # Molar volume of the solute domain
    ov.domain_molar_volume = R.volume.domain * units.Angs3tocm3permol

    # Density of the solute and of the solvent 
    ov.density.solute = R.density.solute * units.SitesperAngs3tomolperL
    ov.density.solvent = R.density.solvent * units.SitesperAngs3tomolperL
    ov.density.solvent_bulk = R.density.solvent_bulk * units.SitesperAngs3tomolperL

    # Solvent molar volume ov.solvent_molar_volume = 1000 / ov.density.solvent
    ov.solvent_molar_volume = 1000 / ov.density.solvent
    ov.solvent_molar_volume_bulk = 1000 / ov.density.solvent_bulk

    # Solute molar volume computed from solvent density in bulk
    if R.autocorrelation
        ov.solute_molar_volume = ov.solvent_molar_volume
    else
        ov.solute_molar_volume =
            units.Angs3tocm3permol *
            (R.density.solvent_bulk * R.volume.total - R.solvent.nmols) /
            R.density.solvent_bulk
    end

    return ov
end
