#=

Structure that contains the detailed input options.

=#
struct Options

    # Frames to be considered
    firstframe::Int
    lastframe::Int
    stride::Int

    # Random sampling options
    irefatom::Int
    n_random_samples::Int

    # Histogram and bulk region options
    binstep::Float64
    cutoff::Float64

    # Linked cell length will be (cutoff/lcell), might be tuned for maximum
    # performance
    lcell::Int

    # Force garbage collection in parallel runs to avoid memory overflow, 
    # whenever free memory in the system is smaller than GC_threshold
    gc::Bool
    gc_threshold::Float64

    # Random number generator options
    seed::Int
    stable_rng::Bool

    # Manually chose on how many threads to run (0 to use all threads)
    nthreads::Int

    # Do not show any output on the screen on execution of mddf routines
    silent::Bool

end

"""
    Options(;
        firstframe::Int = 1,
        lastframe::Int = -1,
        stride::Int = 1,
        irefatom::Int = -1,
        n_random_samples::Union{Nothing,Int} = nothing,
        binstep::Float64 = 0.02,
        cutoff::Union{Nothing,Real} = nothing,
        lcell::Int = 1,
        gc::Bool = true,
        gc_threshold::Float64 = 0.3,
        seed::Int = 321,
        stable_rng::Bool = false,
        nthreads::Int = 0,
        silent::Bool = false
    )

Create an Options object with the specified options. 

`n_random_samples` is the number of random (ideal-gas) samples computed per frame, used to normalize
the distributions. For self-correlations (solute and solvent are the same molecules), each sample is
a configuration in which all molecules are placed with random positions and orientations, and the 
distances of every molecule to all other molecules are computed, as for the actual frame. For 
cross-correlations with more than one solute molecule, each sample is a configuration in which all
solvent molecules are placed with random positions and orientations, and the distances of the solvent
molecules to every solute molecule are computed. For a single solute molecule, each sample is a random 
configuration of the solvent molecules around it. By default (`n_random_samples = nothing`), 1 sample 
per frame is used for self-correlations and for more than one solute molecule, and 10 samples per frame
for a single solute molecule.

The `cutoff` is the maximum distance to the solute for which the distribution functions are 
computed. The bulk density of the solvent is estimated from the solvent molecules beyond the cutoff.
The default value is `cutoff = 10.0`, but it is recommended to set it according to the size of 
the system and to the correlation length of the distribution function.
"""
function Options(;
    firstframe::Int=1,
    lastframe::Int=-1,
    stride::Int=1,
    irefatom::Int=-1,
    n_random_samples::Union{Nothing,Int}=nothing,
    binstep::Float64=0.02,
    cutoff::Union{Nothing,Real}=nothing,
    lcell::Int=1,
    gc::Bool=true,
    gc_threshold::Float64=0.3,
    seed::Int=321,
    stable_rng::Bool=false,
    nthreads::Int=0,
    silent::Bool=false,
    kargs...
)

    _check_invalid_options(kargs)

    # Check for simple input errors
    if stride < 1
        throw(ArgumentError("in MDDF options: stride cannot be less than 1. "))
    end
    if lastframe > 0 && lastframe < firstframe
        throw(ArgumentError("in MDDF options: lastframe must be greater or equal to firstframe. "))
    end
    if !isnothing(n_random_samples) && n_random_samples < 1
        throw(ArgumentError("""\n 
            in MDDF options: n_random_samples must be greater than 0. 
            To skip the normalization of the distribution, use the coordination_number(...) function instead of mddf(...).

        """))
    end

    if isnothing(cutoff)
        silent || @warn """\n
            Using the default value of the cutoff: cutoff = 10.0. 

            It is recommended to set the cutoff according to the size of the system and
            to the correlation length of the distribution function, for example: 
            
                Options(cutoff = 12.0)

        """ _file = nothing _line = nothing
        cutoff = 10.0
    end
    if cutoff <= 0
        throw(ArgumentError("in MDDF options: cutoff must be positive."))
    end
    if (cutoff / binstep) % 1 > 1.e-5
        throw(ArgumentError("in MDDF options: cutoff must be a multiple of binstep."))
    end

    return Options(
        firstframe,
        lastframe,
        stride,
        irefatom,
        # 0 means automatic: the value is set when the computation is initialized (see _set_n_random_samples)
        isnothing(n_random_samples) ? 0 : n_random_samples,
        binstep,
        cutoff,
        lcell,
        gc,
        gc_threshold,
        seed,
        stable_rng,
        nthreads,
        silent
    )
end

#=
    _set_n_random_samples(options::Options, all_to_all::Bool)

Returns a copy of the options with the default number of random samples per frame,
if it was not set by the user: 1 for all-to-all random configurations (self-correlations,
or more than one solute molecule), and 10 for a single solute molecule.

=#
function _set_n_random_samples(o::Options, all_to_all::Bool)
    o.n_random_samples > 0 && return o
    return Options(
        o.firstframe, o.lastframe, o.stride, o.irefatom,
        all_to_all ? 1 : 10,
        o.binstep, o.cutoff, o.lcell, o.gc, o.gc_threshold,
        o.seed, o.stable_rng, o.nthreads, o.silent,
    )
end

#=
    _check_invalid_options(kargs)

Throws informative errors for options removed or renamed in version 4.0.0, and for unknown options.

=#
const _renamed_options = Dict(:GC => :gc, :GC_threshold => :gc_threshold, :StableRNG => :stable_rng)
function _check_invalid_options(kargs)
    for key in keys(kargs)
        if key in (:dbulk, :usecutoff, :bulk_range)
            throw(ArgumentError("""\n
                The `$key` option was removed in version 4.0.0. The bulk density of the solvent 
                is estimated from the region beyond the cutoff, which is set with:

                    Options(cutoff = ...)

            """))
        elseif haskey(_renamed_options, key)
            throw(ArgumentError("""\n
                The `$key` option was renamed to `$(_renamed_options[key])` in version 4.0.0.

            """))
        else
            throw(ArgumentError("""\n
                Unknown option: `$key`. The valid options are:

                    $(join(fieldnames(Options), ", "))

                Type `? Options` for the documentation.

            """))
        end
    end
    return nothing
end

@testitem "Options" begin
    using ComplexMixtures
    o = Options(silent=true)
    @test o.cutoff == 10.0
    o = Options(cutoff=12.0)
    @test o.cutoff == 12.0
    @test_logs (:warn,) Options()
    @test_logs Options(cutoff=12.0)

    # Default number of random samples: 1 for self-correlations, 10 for cross-correlations
    o = Options(cutoff=12.0)
    @test o.n_random_samples == 0
    @test ComplexMixtures._set_n_random_samples(o, true).n_random_samples == 1
    @test ComplexMixtures._set_n_random_samples(o, false).n_random_samples == 10
    o = Options(cutoff=12.0, n_random_samples=5)
    @test ComplexMixtures._set_n_random_samples(o, true).n_random_samples == 5
    @test ComplexMixtures._set_n_random_samples(o, false).n_random_samples == 5

    # Removed, renamed, and unknown options
    @test_throws "removed in version 4.0.0" Options(bulk_range=(10.0, 14.0))
    @test_throws "removed in version 4.0.0" Options(dbulk=10.0)
    @test_throws "removed in version 4.0.0" Options(usecutoff=true)
    @test_throws "renamed to `gc`" Options(GC=true)
    @test_throws "renamed to `gc_threshold`" Options(GC_threshold=0.5)
    @test_throws "renamed to `stable_rng`" Options(StableRNG=true)
    @test_throws "Unknown option: `cutof`" Options(cutof=10.0)

    # input errors
    @test_throws ArgumentError Options(cutoff=10.0, binstep=0.3)
    @test_throws ArgumentError Options(cutoff=-1.0)
    @test_throws ArgumentError Options(stride=0)
    @test_throws ArgumentError Options(lastframe=1, firstframe=2)
    @test_throws ArgumentError Options(n_random_samples=0)
end

function Base.show(io::IO, o::Options)
    print(io,chomp("""
    $bars
    Options - ComplexMixtures 
    $bars
    Trajectory frames: 
        First frame to be considered: firstframe = $(o.firstframe) 
        Last frame to be considered (-1 is last): lastframe = $(o.lastframe)
        Stride: stride = $(o.stride)

    Cutoff and histogram:
        Bin step of histogram: binstep = $(o.binstep)
        Cutoff: cutoff = $(o.cutoff) (bulk density estimated beyond the cutoff)
   
    Computation details: 
        Reference atom for random rotations: irefatom = $(o.irefatom)
        Number of random samples per frame: n_random_samples = $(o.n_random_samples == 0 ? "automatic (1 for self-correlations or more than one solute molecule, 10 otherwise)" : o.n_random_samples)
        Linked cell partition: lcell = $(o.lcell)
        Force garbage collection: gc = $(o.gc)
        Memory threshold for garbage collection: gc_threshold = $(o.gc_threshold)
        Seed for random number generator: seed = $(o.seed)
        Use stable random number generator: stable_rng = $(o.stable_rng) 
        Number of threads to use (0 is all): nthreads = $(o.nthreads)
        Silent output: $(o.silent)
    $bars
    """))
end

#=

Structure that contains data that is trajectory file-specific, to be used within
the Results data structure, particularly when the results are merged from multiple files.

=#
@kwdef mutable struct TrajectoryFileOptions
    const filename::String
    const options::Options
    const irefatom::Int
    const lastframe_read::Int
    nframes_read::Int
    # Statistical weights of each frame of each file in the trajectory. An empty vector
    # means that all frames have the same statistical weight
    const frame_weights::Vector{Float64}
end

function Base.show(io::IO, o::TrajectoryFileOptions)
    bars82 = repeat("-", 82)
    print(io,chomp("""
    $bars82
    Trajectory file: 
        $(o.filename)

    $(o.options)
    $bars82
    Reference atom for random rotations: irefatom = $(o.irefatom)
    Last frame read: lastframe_read = $(o.lastframe_read)
    Number of frames read: nframes_read = $(o.nframes_read)
    Statistical weights of each frame: frame_weights = $(print_vector_summary(o.frame_weights))
    $bars82
    """))
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", vec::Vector{TrajectoryFileOptions})
    print(io, chomp("""
        Vector{ComplexMixtures.TrajectoryFileOptions}
            Number of files: $(length(vec))
            Total number of frames: $(sum(o.nframes_read for o in vec))
        """))
end