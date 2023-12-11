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

    # Statistical weights of each frame in the trajectory. An empty vector
    # means that all frames have the same statistical weight
    frame_weights::Vector{Float64} = Float64[]

    # Groups: a list of groups of atoms for which the contributions will
    # be recorded. Each group is a list of atom indexes. If the list of 
    # groups is empty, the contributions will be recorded for all (type of) atoms   
    solute_groups::Vector{Vector{Int}} = Vector{Int}[]
    solvent_groups::Vector{Vector{Int}} = Vector{Int}[]

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
   
    Groups to compute contributions are pre-defined: 
        For solute: $(isempty(o.solute_groups) ? "no" : "$(length(o.solute_groups)) group(s)")
        For solvent: $(isempty(o.solvent_groups) ? "no" : "$(length(o.solvent_groups)) group(s)")

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
        frame_weights = copy(o.frame_weights),
        solute_groups = Vector{Int}[ copy(g) for g in o.solute_groups ],
        solvent_groups = Vector{Int}[ copy(g) for g in o.solvent_groups ],
    )
end

#
# Merge function for options
#
function Base.merge(O::Vector{Options})
    local options
    empty_solute_groups = false
    empty_solvent_groups = false
    empty_frame_weights = false
    # Check if frame weights were provided to none, all, or some trajectories
    # If there are weights for some trajectories, but not all, issue a warning
    # and fill the empty weights with 1.0
    nframeweights = count(isempty, o.frame_weights for o in O)
    if nframeweights != 0 && nframeweights != length(O)
        @warn begin 
            """
            Frame weights provided only for some results. 
            The frame weights will be empty and should be provided manually or not used for further analysis.
            """ 
        end _file=nothing _line=nothing
        empty_frame_weights = true
    end
    for i in eachindex(O)
        if i == firstindex(O) 
            options = copy(O[i])
            continue
        end
        if !empty_frame_weights
            append!(options.frame_weights, O[i].frame_weights)
        end
        # The group definitions must be the same for all results
        if !empty_solute_groups && !isequal(options.solute_groups, O[i].solute_groups)
            @warn begin
                """
                Groups for solute are not the same for all results. 
                The merged solute_groups will be empty and won't be meaningful.
                """
            end _file=nothing _line=nothing
            empty_solute_groups = true
        end
        if !empty_solvent_groups && !isequal(options.solvent_groups, O[i].solvent_groups)
            @warn begin
                """
                Groups for solvent are not the same for all results. 
                The merged solvent_groups will be empty and won't be meaningful.
                """
            end _file=nothing _line=nothing
            empty_solvent_groups = true
        end
    end
    empty_frame_weights && (empty!(options.frame_weights))
    empty_solute_groups && (empty!(options.solute_groups))
    empty_solvent_groups && (empty!(options.solvent_groups))
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
    @test @test_logs (:warn,) merge([o1, o2]).frame_weights == Float64[]
    o1 = Options(frame_weights=[0.5, 1.0])
    om = merge([o1, o2])
    @test om.frame_weights == [0.5, 1.0, 1.0, 2.0] 

    o1 = Options(solute_groups=[[1,2,3]])
    @test o1.solute_groups == [[1,2,3]]
    @test o1.solvent_groups == Vector{Int}[]
    o2 = Options(solute_groups=[[1,2,3]])
    om = merge([o1, o2])
    @test om.solute_groups == [[1,2,3]]
    @test om.solvent_groups == Vector{Int}[]
    o2 = Options(solute_groups=[[1,2,3], [4,5,6]])
    om = @test_logs (:warn, ) merge([o1, o2])
    @test om.solute_groups == Vector{Int}[]
    @test om.solvent_groups == Vector{Int}[]
    o1 = Options(solute_groups=[[1,2,3], [4,5,6]])
    o2 = Options(solute_groups=[[1,2,3]])
    om = @test_logs (:warn, ) merge([o1, o2])
    @test om.solute_groups == Vector{Int}[]
    @test om.solvent_groups == Vector{Int}[]

    o1 = Options(solvent_groups=[[1,2,3]])
    @test o1.solvent_groups == [[1,2,3]]
    @test o1.solute_groups == Vector{Int}[]
    o2 = Options(solvent_groups=[[1,2,3]])
    om = merge([o1, o2])
    @test om.solvent_groups == [[1,2,3]]
    @test om.solute_groups == Vector{Int}[]
    o1 = Options(solvent_groups=[[1,2,3]])
    o2 = Options(solvent_groups=[[1,2,3], [4,5,6]])
    om = @test_logs (:warn, ) merge([o1, o2])
    @test om.solvent_groups == Vector{Int}[]
    @test om.solute_groups == Vector{Int}[]
    o1 = Options(solvent_groups=[[1,2,3], [4,5,6]])
    o2 = Options(solvent_groups=[[1,2,3]])
    om = @test_logs (:warn, ) merge([o1, o2])
    @test om.solvent_groups == Vector{Int}[]
    @test om.solute_groups == Vector{Int}[]

    o1 = Options(solvent_groups=[[1,2,3]], solute_groups=[[4,5,6]])
    o2 = Options(solvent_groups=[[1,2,3]], solute_groups=[[4,5,6]])
    om = merge([o1, o2])
    @test om.solvent_groups == [[1,2,3]]
    @test om.solute_groups == [[4,5,6]]

end
