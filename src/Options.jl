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
    dbulk::Float64
    cutoff::Float64
    usecutoff::Bool

    # Linked cell length will be (cutoff/lcell), might be tunned for maximum
    # performance
    lcell::Int

    # Force garbage collection in parallel runs to avoid memory overflow, 
    # whenever free memory in the system is smaller than GC_threshold
    GC::Bool
    GC_threshold::Float64

    # Random number generator options
    seed::Int
    StableRNG::Bool

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
        n_random_samples::Int = 10,
        binstep::Float64 = 0.02,
        dbulk::Union{Nothing,Real} = nothing,
        cutoff::Union{Nothing,Real} = nothing,
        usecutoff::Union{Nothing,Bool} = nothing,
        bulk_range=nothing,
        lcell::Int = 1,
        GC::Bool = true,
        GC_threshold::Float64 = 0.3,
        seed::Int = 321,
        StableRNG::Bool = false,
        nthreads::Int = 0,
        silent::Bool = false
    )

Create an Options object with the specified options. 

"""
function Options(;
    firstframe::Int = 1,
    lastframe::Int = -1,
    stride::Int = 1,
    irefatom::Int = -1,
    n_random_samples::Int = 10,
    binstep::Float64 = 0.02,
    dbulk::Union{Nothing,Real} = nothing,
    cutoff::Union{Nothing,Real} = nothing,
    usecutoff::Union{Nothing,Bool} = nothing,
    bulk_range=nothing,
    lcell::Int = 1,
    GC::Bool = true,
    GC_threshold::Float64 = 0.3,
    seed::Int = 321,
    StableRNG::Bool = false,
    nthreads::Int = 0,
    silent::Bool = false
)

    # warning flag for default values of dbulk, cutoff, and usecutoff
    warn = false

    # Check for simple input errors
    if stride < 1
       throw(ArgumentError("in MDDF options: stride cannot be less than 1. "))
    end
    if lastframe > 0 && lastframe < firstframe
        throw(ArgumentError("in MDDF options: lastframe must be greater or equal to firstframe. "))
    end


    if !isnothing(bulk_range) && any(!isnothing, (dbulk, cutoff, usecutoff))
        throw(ArgumentError("""\n
            The bulk_range argument implies that dbulk, cutoff, and usecutoff are not needed. 
        
        """)) 
    end
    if all(isnothing, (bulk_range, dbulk, cutoff, usecutoff))
        dbulk = 10.0
        cutoff = 10.0
        usecutoff = false
        warn = true
    elseif !isnothing(bulk_range)
        if length(bulk_range) != 2
            throw(ArgumentError("""\n
                bulk_range must be a tuple or vector with two elements, corresponding to dbulk and cutoff.
                Example: Options(;bulk_range = (8.0, 12.0))

            """))
        end
        dbulk, cutoff = bulk_range
        usecutoff = true
    else
        if isnothing(dbulk)
            dbulk = 10.0
            warn = true
        end
        if isnothing(usecutoff)
            usecutoff = false
            warn = true
        end
        if isnothing(cutoff)
            if usecutoff
                cutoff = dbulk + 4.0
                warn = true
            else
                warn = true
                cutoff = dbulk
            end
        else
            if !usecutoff
                throw(ArgumentError("in MDDF options: cutoff was defined with usecutoff set to false"))
            end
        end
    end
    if warn && !silent
        @warn """\n
            Using default values for dbulk, cutoff and/or usecutoff: 
            
                dbulk = $(dbulk)
                cutoff = $(cutoff)
                usecutoff = $(usecutoff)
            
            It is recommended to set bulk_range manually, according to the system size and correlations of the distribution function.

        """ _file=nothing _line=nothing
    end
    if usecutoff && dbulk >= cutoff
        throw(ArgumentError(" in MDDF options: The bulk volume is zero (dbulk must be smaller than cutoff). "))
    end
    if (cutoff / binstep) % 1 > 1.e-5
        throw(ArgumentError("in MDDF options: cutoff must be a multiple of binstep."))
    end
    if (dbulk / binstep) % 1 > 1.e-5
        throw(ArgumentError("in MDDF options: dbulk must be a multiple of binstep."))
    end

    return Options(
        firstframe,
        lastframe,
        stride,
        irefatom,
        n_random_samples,
        binstep,
        dbulk,
        cutoff,
        usecutoff,
        lcell,
        GC,
        GC_threshold,
        seed,
        StableRNG,
        nthreads,
        silent
    )
end

@testitem "Options" begin
    using ComplexMixtures
    o = Options()
    @test o.dbulk == 10.0
    @test o.cutoff == 10.0
    @test o.usecutoff == false
    o = Options(bulk_range = (10.0, 14.0))
    @test o.dbulk == 10.0
    @test o.cutoff == 14.0
    @test o.usecutoff == true
    o = Options(dbulk=10.0)
    @test o.dbulk == 10.0
    @test o.cutoff == 10.0
    @test o.usecutoff == false 
    o = Options(dbulk = 10.0, usecutoff = false)
    @test o.dbulk == 10.0
    @test o.cutoff == 10.0
    @test o.usecutoff == false
    o = Options(bulk_range = (10.0, 14.0))
    @test o.dbulk == 10.0
    @test o.cutoff == 14.0
    @test o.usecutoff == true

    # input errors
    @test_throws ArgumentError Options(dbulk = 10.0, binstep = 0.3)
    @test_throws ArgumentError Options(dbulk = 10.0, cutoff = 10.0, usecutoff = true)
    @test_throws ArgumentError Options(dbulk = 6.0, cutoff = 10.0, binstep = 0.3, usecutoff = true)
    @test_throws ArgumentError Options(dbulk = 8.0, cutoff = 10.0, binstep = 0.3, usecutoff = true)
    @test_throws ArgumentError Options(bulk_range = (10.0, 12.0), binstep = 0.3)
    @test_throws ArgumentError Options(stride = 0)
    @test_throws ArgumentError Options(lastframe = 1, firstframe = 2)
    @test_throws ArgumentError Options(bulk_range = (12.0, 10.0))
    @test_throws ArgumentError Options(bulk_range = (12.0, 10.0), dbulk=8.0)
    @test_throws ArgumentError Options(bulk_range = (12.0, 10.0), cutoff=15.0)
    @test_throws ArgumentError Options(bulk_range = (12.0, 10.0), usecutoff=false)
    @test_throws ArgumentError Options(dbulk = 10.0, cutoff = 15.0, usecutoff = false)

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
        Bulk range: $(o.usecutoff ? "$(o.dbulk) - $(o.cutoff)" : ">= $(o.dbulk)")
        (dbulk = $(o.dbulk), cutoff = $(o.cutoff), usecutoff = $(o.usecutoff))
   
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