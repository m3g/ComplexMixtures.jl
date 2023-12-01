"""

$(TYPEDEF)

Structure that contains the detailed input options.

$(TYPEDFIELDS)

"""
Base.@kwdef struct Options

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

    # Statistical weights of each frame in the trajectory. An empty vector
    # means that all frames have the same statistical weight
    frame_weights::Vector{Float64} = Float64[]

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

    Statistical weight of frames:
        frame_weights = $(
            isempty(o.frame_weights) ? 1.0 : 
                if length(o.frame_weights) > 4
                    "["*join(round.(o.frame_weights[begin:begin+1],digits=2), ", ")*
                    ", ... ,"*
                    join(round.(o.frame_weights[end-1:end], digits=2), ", ")*" ]"
                else
                   o.frame_weights 
                end
            )
        length of weights vector = $(length(o.frame_weights))

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
    """
    ))
    return nothing
end

#
# We need to define a custom isequal function because now (as of 1.4.0)
# the Options struct has a mutable field, which makes the comparison 
# of the complete structs evaluate to false
#
import Base: == 
function ==(o1::Options, o2::Options)
    eq = true
    for field in fieldnames(Options)
        eq = getfield(o1, field) == getfield(o2, field)
    end
    return eq
end
