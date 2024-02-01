"""

$(TYPEDEF)

Structure that contains the detailed input options.

$(TYPEDFIELDS)

"""
@kwdef struct Options

    firstframe::Int = 1
    lastframe::Int = -1
    stride::Int = 1

    irefatom::Int = -1
    n_random_samples::Int = 10

    binstep::Float64 = 0.02
    dbulk::Float64 = 10.0
    cutoff::Float64 = 10.0
    usecutoff::Bool = false

    # Linked cell length will be (cutoff/lcell), might be tunned for maximum
    # performance
    lcell::Int = 1

    # Force garbage collection in parallel runs to avoid memory overflow, 
    # whenever free memory in the system is smaller than GC_threshold
    GC::Bool = true
    GC_threshold::Float64 = 0.3

    # Random number generator options
    seed::Int = 321
    StableRNG::Bool = false

    # Manually chose on how many threads to run (0 to use all threads)
    nthreads::Int = 0

    # Do not show any output on the screen on execution of mddf routines
    silent::Bool = false

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

    Bulk region, cutoff, and histogram:
        Bin step of histogram: binstep = $(o.binstep)
        Bulk distance: dbulk = $(o.dbulk)
        Use cutoff diffrent from dbulk: usecuoff = $(o.usecutoff)
        Cutoff: cutoff = $(o.cutoff)
   
    Computation details: 
        Reference atom for random rotations: irefatom = $(o.irefatom)
        Number of random samples per frame: n_random_samples = $(o.n_random_samples)
        Linked cell partition: lcell = $(o.lcell)
        Force garbage collection: GC = $(o.GC)
        Memory threshold for GC: GC_threshold = $(o.GC_threshold)
        Seed for random number generator: $(o.seed)
        Use stable random number generator: StableRNG = $(o.StableRNG) 
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