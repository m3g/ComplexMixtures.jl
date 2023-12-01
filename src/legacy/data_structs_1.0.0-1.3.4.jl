#
# Legacy data structures, from 1.0.0 to 1.3.4
#
Base.@kwdef struct LegacyOptionsA

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

struct LegacySolSummaryA
    natoms::Int
    nmols::Int
    natomspermol::Int
end

@kwdef mutable struct LegacyDensityA
    solute::Float64 = 0.0
    solvent::Float64 = 0.0
    solvent_bulk::Float64 = 0.0
end

@kwdef mutable struct LegacyVolumeA
    total::Float64
    bulk::Float64
    domain::Float64
    shell::Vector{Float64}
end

@kwdef mutable struct LegacyResultA{T<:VecOrMat{Float64}}

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
    mddf::Vector{Float64} = zeros(nbins)
    kb::Vector{Float64} = zeros(nbins)

    # Properties of the solute and solvent selections
    autocorrelation::Bool
    solvent::LegacySolSummaryA
    solute::LegacySolSummaryA

    # Atomic contributions to the MDDFs
    solute_atom::T = zeros(nbins, solute.natomspermol)
    solvent_atom::T = zeros(nbins, solvent.natomspermol)

    # Data to compute a RDF and the KB integral from this count
    rdf_count::Vector{Float64} = zeros(nbins)
    rdf_count_random::Vector{Float64} = zeros(nbins)
    sum_rdf_count::Vector{Float64} = zeros(nbins)
    sum_rdf_count_random::Vector{Float64} = zeros(nbins)
    rdf::Vector{Float64} = zeros(nbins)
    kb_rdf::Vector{Float64} = zeros(nbins)

    # Overall densities and volumes
    density::LegacyDensityA = LegacyDensityA()
    volume::LegacyVolumeA = LegacyVolumeA(nbins)

    # Options of the calculation
    options::LegacyOptionsA
    irefatom::Int
    lastframe_read::Int
    nframes_read::Int

    # File name(s) of the trajectories in this results 
    files::Vector{String}
    weights::Vector{Float64}
end

#
# How to update a result read with this data structure to the current one
#
function update_result_structure(r::LegacyResultA)
    options = Options(
        firstframe = r.options.firstframe,
        lastframe = r.options.lastframe,
        stride = r.options.stride,
        irefatom = r.options.irefatom,
        n_random_samples = r.options.n_random_samples,
        binstep = r.options.binstep,
        dbulk = r.options.dbulk,
        cutoff = r.options.cutoff,
        usecutoff = r.options.usecutoff,
        lcell = r.options.lcell,
        GC = r.options.GC,
        GC_threshold = r.options.GC_threshold,
        seed = r.options.seed,
        StableRNG = r.options.StableRNG,
        nthreads = r.options.nthreads,
        silent = r.options.silent,
        frame_weights = Float64[], # new field in 1.3.5
    )
    solvent = SolSummary(r.solvent.natoms, r.solvent.nmols, r.solvent.natomspermol)
    solute = SolSummary(r.solute.natoms, r.solute.nmols, r.solute.natomspermol)
    density = Density(r.density.solute, r.density.solvent, r.density.solvent_bulk)
    volume = Volume(r.volume.total, r.volume.bulk, r.volume.domain, r.volume.shell)
    result = Result(
        Version = pkgversion(@__MODULE__), # new field in 1.3.5
        nbins = r.nbins,
        dbulk = r.dbulk,
        cutoff = r.cutoff,
        d = r.d,
        md_count = r.md_count,
        md_count_random = r.md_count_random,
        coordination_number = r.coordination_number, 
        coordination_number_random = r.coordination_number_random,
        mddf = r.mddf,
        kb = r.kb,
        autocorrelation = r.autocorrelation,
        solvent = solvent,
        solute = solute,
        solute_atom = r.solute_atom,
        solvent_atom = r.solvent_atom,
        rdf_count = r.rdf_count,
        rdf_count_random = r.rdf_count_random,
        sum_rdf_count = r.sum_rdf_count,
        sum_rdf_count_random = r.sum_rdf_count_random,
        rdf = r.rdf,
        kb_rdf = r.kb_rdf,
        density = density,
        volume = volume,
        options = options,
        irefatom = r.irefatom,
        lastframe_read = r.lastframe_read,
        nframes_read = r.nframes_read,
        files = r.files,
        weights = r.weights,
    )
    return result
end

#
# Functions to save the results to a file
#
StructTypes.StructType(::Type{LegacySolSummaryA}) = StructTypes.Struct()
StructTypes.StructType(::Type{LegacyResultA{Vector{Float64}}}) = StructTypes.Struct()
StructTypes.StructType(::Type{LegacyResultA{Matrix{Float64}}}) = StructTypes.Struct()
StructTypes.StructType(::Type{LegacyDensityA}) = StructTypes.Struct()
StructTypes.StructType(::Type{LegacyVolumeA}) = StructTypes.Struct()
StructTypes.StructType(::Type{LegacyOptionsA}) = StructTypes.Struct()

# Read legacy result json files
function load_json_LegacyA(filename)
    f = open(filename, "r")
    R = JSON3.read(f, LegacyResultA{Vector{Float64}})
    close(f)
    # Need to reshape the solute and solvent atom contributions, because the data is read in a single column
    solute_atom = reshape(R.solute_atom, R.nbins, :)
    solvent_atom = reshape(R.solvent_atom, R.nbins, :)
    r_names = fieldnames(LegacyResultA)
    r = LegacyResultA{Matrix{Float64}}(
        ntuple(length(r_names)) do i
            r_names[i] == :solute_atom ? solute_atom :
            r_names[i] == :solvent_atom ? solvent_atom : getfield(R, r_names[i])
        end...,
    )
    return update_result_structure(r)
end
