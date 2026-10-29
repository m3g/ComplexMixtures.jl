#=

$(TYPEDEF)

**LEGACY** support for versions 2.18.0 - 3.1.0. This data structure will read 
the result of a legacy file and be used to construct an updated version of the 
result data structure. In v4.0.0 the data structure was renamed from `Result` to
`MDDFResult`, and the `rdf*` fields were renamed to `proximal*`.

$(TYPEDFIELDS)

=#
@kwdef mutable struct Result_3_1_0

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
    # of using JSON, we use Vector{Vector{Float64}}, which is 
    # read directly.
    solute_group_count::Vector{Vector{Float64}}
    solvent_group_count::Vector{Vector{Float64}}
    solute_group_count_random::Vector{Vector{Float64}}
    solvent_group_count_random::Vector{Vector{Float64}}

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

function _load(filename, ::Type{Result_3_1_0})
    R_old = try
        open(filename, "r") do io
            JSON.parse(io, Result_3_1_0)
        end
    catch
        throw(ArgumentError("""\n 
            The file $filename does not appear to contain a MDDFResult object.

        """))
    end
    # Convert legacy Result_3_1_0 to current MDDFResult data structure
    # This conversion must be updated if a new MDDFResult data structure is set.
    return MDDFResult(;
        # The converted data structure is in the format of the current version
        Version = pkgversion(@__MODULE__),
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
        # Group contributions
        solute_group_count = R_old.solute_group_count,
        solvent_group_count = R_old.solvent_group_count,
        solute_group_count_random = R_old.solute_group_count_random,
        solvent_group_count_random = R_old.solvent_group_count_random,
        # Proximal distribution (`rdf*` fields in versions < 4.0.0)
        proximal_count = R_old.rdf_count,
        proximal_count_random = R_old.rdf_count_random,
        sum_proximal_count = R_old.sum_rdf_count,
        sum_proximal_count_random = R_old.sum_rdf_count_random,
        proximal = R_old.rdf,
        kb_proximal = R_old.kb_rdf,
        # Overall densities and volumes
        density = R_old.density,
        volume = R_old.volume,
        # Trajectory files and weights
        files = R_old.files,
        weights = R_old.weights,
    )
end

@testitem "Result_3_1_0 reading" begin
    using ComplexMixtures
    using LazyArtifacts
    using ComplexMixtures: test_file
    file = test_file("glycerol_protein_glyc.json")
    @test ComplexMixtures._get_version(file) < v"4.0.0-DEV"
    R = load(MDDFResult, file)
    raw = ComplexMixtures.JSON.parsefile(file)
    @test R.Version == pkgversion(ComplexMixtures)
    @test R.mddf ≈ Float64.(raw["mddf"])
    @test R.proximal ≈ Float64.(raw["rdf"])
    @test R.proximal_count ≈ Float64.(raw["rdf_count"])
    @test R.proximal_count_random ≈ Float64.(raw["rdf_count_random"])
    @test R.sum_proximal_count ≈ Float64.(raw["sum_rdf_count"])
    @test R.sum_proximal_count_random ≈ Float64.(raw["sum_rdf_count_random"])
    @test R.kb_proximal ≈ Float64.(raw["kb_rdf"])
    @test length(R.solute_group_count_random) == length(raw["solute_group_count_random"])
    # The converted result is saved in the current format, and can be reloaded
    tmp = tempname()
    save(tmp, R)
    @test !occursin("\"rdf\"", read(tmp, String))
    @test load(MDDFResult, tmp) == R
    rm(tmp)
end
