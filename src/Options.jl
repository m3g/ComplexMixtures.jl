"""

$(TYPEDEF)

Structure that contains the detailed input options.

$(TYPEDFIELDS)

"""
@with_kw_noshow struct Options

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

# A copy function is needed for the merge function. Could be deepcopy,
# but feels safer to control the copy here: need to copy the mutable
# fields explicitly
import Base: copy
function copy(o::Options)
    return Options(
        firstframe = o.firstframe,
        lastframe = o.lastframe,
        stride = o.stride,
        n_random_samples = o.n_random_samples,
        binstep = o.binstep,
        dbulk = o.dbulk,
        cutoff = o.cutoff,
        usecutoff = o.usecutoff,
        lcell = o.lcell,
        GC = o.GC,
        GC_threshold = o.GC_threshold,
        seed = o.seed,
        StableRNG = o.StableRNG,
        nthreads = o.nthreads,
        silent = o.silent,
        frame_weights = copy(o.frame_weights)
    )
end

#
# Merge function for options
#
function Base.merge(O::Vector{Options})
    options = copy(O[1])
    for i in 2:length(O)
        if (isempty(options.frame_weights) && !isempty(O[i].frame_weights)) ||
           (!isempty(options.frame_weights) && isempty(O[i].frame_weights))
            throw(ArgumentError("Frame weights provided only for some results. The merged frame_weights will be empty and won't be meaningful."))
        else
            append!(options.frame_weights, O[i].frame_weights)
        end
    end
    return options
end

@testitem "copy/merge Options" begin
    o1 = Options()
    o2 = copy(o1)
    @test o1 == o2
    o2 = Options()
    om = merge([o1, o2])
    @test om == o1
    o2 = Options(frame_weights=[1.0, 2.0])
    @test_throws ArgumentError merge([o1, o2])
    o1 = Options(frame_weights=[0.5, 1.0])
    om = merge([o1, o2])
    @test om.frame_weights == [0.5, 1.0, 1.0, 2.0] 
end
