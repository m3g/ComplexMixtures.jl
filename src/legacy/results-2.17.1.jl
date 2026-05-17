#=

$(TYPEDEF)

**LEGACY** support for versions 2.0.0 - 2.17.1. This data structure will read 
the result of a legacy file and be used to construct an updated version of the 
result data structure.

Structure to contain the results of the MDDF calculation.

$(TYPEDFIELDS)

The Result{Vector{Float64}} parametric type is necessary only for reading the JSON3 saved file. 

=#
@kwdef mutable struct Result_2_17_1

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

function load(filename, ::Type{Result_2_17_1})
    R_old = try
        open(filename, "r") do io
            JSON3.read(io, Result_2_17_1)
        end
    catch
        throw(ArgumentError("""\n 
            The file $filename does not appear to contain a Result object.

        """))
    end
    # Convert legacy Result_2.17_1 to current Result data structure
    # This conversion must be updated if a new Result data structure is set.
    return Result(;
        # ComplexMixtures version that generated this results
        Version = R_old.Version,
        # Histogram properties
        nbins = R_old.nbins,
        dbulk = R_old.dbulk,
        cutoff = R_old.cutoff,
        d = R_old.d,
        # Data to compute the MDDF distribution and corresponding KB integral
        md_count = R_old.md_count,
        md_count_random = R_old.md_count_random,
        coordination_number = R_old.coordination_number,
        coordination_number_random = R_old.coordination_number_random,
        # The resulting MDDF and kb integrals
        mddf = R_old.mddf,
        kb = R_old.kb,
        # Properties of the solute and solvent selections
        autocorrelation = R_old.autocorrelation,
        solute = R_old.solute, 
        solvent = R_old.solvent,

        solute_group_count = R_old.solute_group_count,
        solvent_group_count = R_old.solvent_group_count,
        solute_group_count_random = Vector{Vector{Float64}}[], # empty in v"2.17.1"
        solvent_group_count_random = Vector{Vector{Float64}}[], # empty in v"2.17.1"

        # Data to compute a RDF and the KB integral from this count
        rdf_count = R_old.rdf_count,
        rdf_count_random = R_old.rdf_count_random,
        sum_rdf_count = R_old.sum_rdf_count,
        sum_rdf_count_random = R_old.sum_rdf_count_random,
        rdf = R_old.rdf,
        kb_rdf = R_old.kb_rdf,

        # Overall densities and volumes
        density = R_old.density,
        volume = R_old.volume,

        #
        # If multiple files, the following vector defines the weight of each file
        # in the final merged results. It is proportional to the number of frames
        # and weights provided for the frames in each file. 
        #
        files = R_old.files,
        weights = R_old.weights,
    )
end

@testitem "Result_2_17_1 reading" begin
    using ComplexMixtures
    using ComplexMixtures: data_dir
    r_new = load(joinpath(data_dir, "NAMD/tmao_tmao.json"))
    r_old = load(joinpath(data_dir, "legacy/tmao_tmao.json"))
    @test r_new.mddf ≈ r_old.mddf
    @test r_old.solute_group_count_random == Vector{Float64}[]
    @test r_old.solvent_group_count_random == Vector{Float64}[]
end
