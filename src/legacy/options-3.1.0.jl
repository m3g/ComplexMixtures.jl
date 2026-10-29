#=

**LEGACY** support for reading the options of the results saved with versions 2.0.0 - 3.1.0. 
In version 4.0.0, the `dbulk` and `usecutoff` fields were removed, and the `GC`, `GC_threshold`, 
and `StableRNG` fields were renamed to `gc`, `gc_threshold`, and `stable_rng`.

=#
struct Options_3_1_0
    firstframe::Int
    lastframe::Int
    stride::Int
    irefatom::Int
    n_random_samples::Int
    binstep::Float64
    dbulk::Float64
    cutoff::Float64
    usecutoff::Bool
    lcell::Int
    GC::Bool
    GC_threshold::Float64
    seed::Int
    StableRNG::Bool
    nthreads::Int
    silent::Bool
end

_convert(o::Options_3_1_0) = Options(
    o.firstframe, o.lastframe, o.stride, o.irefatom, o.n_random_samples,
    o.binstep, o.cutoff, o.lcell, o.GC, o.GC_threshold,
    o.seed, o.StableRNG, o.nthreads, o.silent,
)

@kwdef struct TrajectoryFileOptions_3_1_0
    filename::String
    options::Options_3_1_0
    irefatom::Int
    lastframe_read::Int
    nframes_read::Int
    frame_weights::Vector{Float64}
end

_convert(o::TrajectoryFileOptions_3_1_0) = TrajectoryFileOptions(
    filename=o.filename,
    options=_convert(o.options),
    irefatom=o.irefatom,
    lastframe_read=o.lastframe_read,
    nframes_read=o.nframes_read,
    frame_weights=o.frame_weights,
)
