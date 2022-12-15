"""

$(TYPEDEF)

Structure that contains the detailed input options.

$(TYPEDFIELDS)

"""
@with_kw struct Options

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
    GC_threshold::Float64 = 0.1

    # Random number generator options
    seed::Int = 321
    StableRNG::Bool = false

    # Manually chose on how many threads to run (0 to use all threads)
    nthreads::Int = 0

    # Do not show any output on the screen on execution of mddf routines
    silent::Bool = false

end
