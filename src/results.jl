"""

$(TYPEDEF)

Structure to contain the density values obtained from the calculation.

$(TYPEDFIELDS)

"""
@kwdef mutable struct Density
    solute::Float64 = 0.0
    solvent::Float64 = 0.0
    solvent_bulk::Float64 = 0.0
end

function Base.show(io::IO, d::Density)
    println(io, "Density of solute: ", d.solute)
    println(io, "Density of solvent: ", d.solvent)
    print(io, "Density of solvent in bulk: ", d.solvent_bulk)
end

#=
    setbin(d,step)

Function that sets to which histogram bin a data point pertains simple, but important to keep consistency over all calls.

=#
setbin(d, step) = max(1, ceil(Int, d / step))

#=

Structures to contain the volumes obtained from calculations.

=#
@kwdef mutable struct Volume
    total::Float64
    bulk::Float64
    domain::Float64
    shell::Vector{Float64}
end
Volume(nbins::Int) = Volume(0.0, 0.0, 0.0, zeros(Float64, nbins))

function Base.show(io::IO, v::Volume)
    println(io, "Total volume: ", v.total)
    println(io, "Bulk volume: ", v.bulk)
    println(io, "Domain volume: ", v.domain)
    print(io, "Shell volumes: ", v.shell)
end

"""

$(TYPEDEF)

Structure to contain the results of the MDDF calculation.

$(TYPEDFIELDS)

The Result{Vector{Float64}} parametric type is necessary only for reading the JSON3 saved file. 

"""
@kwdef mutable struct Result

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
    # compute group contributions to the MDDFs
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

#
# Initialize the data structure that is returned from the computation, and checks some
# input parameters for consistency
#
function Result(
    trajectory::Trajectory,
    options::Options; 
    trajectory_data::TrajectoryMetaData = TrajectoryMetaData(trajectory, options),
    frame_weights::Vector{Float64} = Float64[]
)

    # Number of bins of the histogram
    nbins = setbin(options.cutoff, options.binstep)

    # If frame weights are provided, the length of the weights vector has to at least
    # of the of number of the last frame to be read
    if !isempty(frame_weights)
        if length(frame_weights) < trajectory_data.lastframe_read
            throw(ArgumentError(chomp("""\n
            The length of the frame_weights vector provided must at least the number of frames to be read.
            
                Input given: length(frame_weights) = $(length(frame_weights))
                             last frame to be read: $(trajectory_data.lastframe_read)
            
            """)))
        end
    else
        frame_weights = ones(trajectory_data.lastframe_read)
    end

    # Initialize group count arrays
    solute_group_count = [ zeros(nbins) for _ in 1:trajectory_data.n_groups_solute ]
    solvent_group_count = [ zeros(nbins) for _ in 1:trajectory_data.n_groups_solvent ]

    return Result(
        nbins = nbins,
        dbulk = options.dbulk,
        cutoff = options.cutoff,
        autocorrelation = isautocorrelation(trajectory),
        solute = trajectory.solute,
        solvent = trajectory.solvent,
        solute_group_count = solute_group_count,
        solvent_group_count = solvent_group_count, 
        files = [
            TrajectoryFileOptions(
                filename = trajectory.filename,
                options = options,
                irefatom = trajectory_data.irefatom, 
                lastframe_read = trajectory_data.lastframe_read,
                nframes_read = trajectory_data.nframes_read,
                frame_weights = frame_weights,
            )
        ],
        weights = [1.0],
    )
end

@testitem "input: argument errors" begin
    using ComplexMixtures: mddf, Trajectory, Options, AtomSelection
    using ComplexMixtures.Testing: data_dir
    using PDBTools: readPDB, select
    atoms = readPDB("$data_dir/NAMD/structure.pdb")
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol = 14)
    water = AtomSelection(select(atoms, "water"), natomspermol = 3)
    trajectory = Trajectory("$data_dir/NAMD/trajectory.dcd", tmao, water)
    @test_throws ArgumentError mddf(trajectory, Options(lastframe = 100))
    @test_throws ArgumentError mddf(trajectory, Options(irefatom = 1000))
    @test_throws ArgumentError mddf(trajectory, Options(lastframe = 5), frame_weights = [1.0, 1.0, 1.0])
end

#
# What to show at the REPL
#
Base.show(io::IO, R::Result) = show(io, overview(R))

#=
Type of calculation: auto or cross-correlations - have different sampling schemes

Choose self and cross versions

In both the self and (cross) non-self cases, the number of random samples is 
n_random_samples*solvent.nmols. However, in the non-self distribution
the sampling of solvent distances is proportional to the number of solute
molecules, thus the md count has to be averaged by solute.nmols. In the
case of the self-distribution, we compute n(n-1)/2 distances, and we will
divide this by n random samples, which is the sampling of the random
distribution. Therefore, we must weight the self-distance count by dividing
it by (n-1)/2, so that we have a count proportional to n as well, leading
to the correct weight relative to the random sample. 
=#
function set_samples(R::Result)
    if R.autocorrelation
        samples = (solvent_nmols = R.solvent.nmols - 1, random = R.files[1].options.n_random_samples)
    else
        samples = (solvent_nmols = R.solvent.nmols, random = R.files[1].options.n_random_samples)
    end
    return samples
end

# autocorrelation can be obtained from the comparison of solute and solvent indices
isautocorrelation(solute_indices::Vector{Int}, solvent_indices::Vector{Int}) =
    solute_indices == solvent_indices ? true : false
isautocorrelation(trajectory::Trajectory) =
    isautocorrelation(trajectory.solute.indices, trajectory.solvent.indices)

#
# Functions to compute volumes of shells
#
#=
    sphericalshellvolume(i,step)

Computes the volume of the spherical shell defined within [(i-1)*step,i*step].

=#
function sphericalshellvolume(i, step)
    rmin = (i - 1) * step
    return (4 * pi / 3) * ((rmin + step)^3 - rmin^3)
end

@testitem "sphericalshellvolume" begin
    using ComplexMixtures: sphericalshellvolume
    @test sphericalshellvolume(1, 1.0) ≈ 4 * pi / 3
    @test sphericalshellvolume(2, 1.0) ≈ 4 * pi / 3 * (8 - 1)
    @test sphericalshellvolume(3, 1.0) ≈ 4 * pi / 3 * (27 - 8)
end

#=
    shellradius(i,step)

Compute the point in which the radius comprises half of the volume of the shell.

=#
function shellradius(i, step)
    rmin = (i - 1) * step
    return (0.5 * ((rmin + step)^3 + rmin^3))^(1 / 3)
end

@testitem "shellradius" begin
    using ComplexMixtures: shellradius
    @test shellradius(1, 0.1) ≈ 0.07937005259840998
    @test shellradius(5, 0.3) ≈ 1.3664650373440481
end

#
# Function to compute the normalization of the weights, given the optional
# frame weights and the number of frames read, which is used in the case
# that frame weights were not provided
#
function sum_frame_weights(R::Result)
    Q = 0.0
    for file in R.files
        i = file.options.firstframe
        s = file.options.stride
        l = file.lastframe_read
        Q += sum(file.frame_weights[i] for i in i:s:l)
    end
    return Q
end

#=
    finalresults!(R::Result, options::Options; coordination_number_only)

Function that computes the final results of all the data computed by averaging according to the sampling of each type of data, and converts to common units.

Computes also the final distribution functions and KB integrals.

This function modified the values contained in the R data structure.

If `coordination_number_only` is true, do not print a warning if there are zero samples in the ideal-gas histogram bins.

=#
function finalresults!(R::Result, options::Options; coordination_number_only)
    R = if !coordination_number_only
        _mddf_final_results!(R, options)
    else
        _coordination_number_final_results!(R, options)
    end
    return R
end

function _mddf_final_results!(R::Result, options::Options)

    # Setup the distance vector
    R.d .= shellradius.(1:R.nbins, options.binstep)

    # Sampling scheme depending on the type of calculation
    samples = set_samples(R)

    # Normalization of number of frames: sum of weights for all frames read
    Q = sum_frame_weights(R)

    # Scale counters by number of samples and frames
    @. R.md_count = R.md_count / (R.solute.nmols * Q)
    @. R.solute_group_count = R.solute_group_count / (R.solute.nmols * Q)
    if R.autocorrelation
        R.solvent_group_count .= R.solute_group_count
    else
        @. R.solvent_group_count = R.solvent_group_count / (R.solute.nmols * Q)
    end
    @. R.md_count_random = R.md_count_random / (samples.random * Q)
    @. R.rdf_count = R.rdf_count / (R.solute.nmols * Q)
    @. R.rdf_count_random = R.rdf_count_random / (samples.random * Q)

    # Volume of each bin shell and of the solute domain
    R.volume.total = R.volume.total / Q
    @. R.volume.shell = R.volume.total * (R.rdf_count_random / samples.solvent_nmols)

    # Solute domain volume
    ibulk = setbin(R.dbulk + 0.5 * R.files[1].options.binstep, R.files[1].options.binstep)
    R.volume.domain = sum(@view(R.volume.shell[1:ibulk-1]))

    # Bulk volume and density properties: either the bulk is considered everything
    # that is not the domain, or the bulk is the region between d_bulk and cutoff,
    # if R.files[1].options.usecutoff is true (meaning that there is a cutoff different from
    # that of the bulk distance)
    if !R.files[1].options.usecutoff
        R.volume.bulk = R.volume.total - R.volume.domain
        n_solvent_in_bulk = samples.solvent_nmols - sum(R.rdf_count)
    else
        n_solvent_in_bulk = sum(@view(R.rdf_count[ibulk:R.nbins]))
        R.volume.bulk = sum(@view(R.volume.shell[ibulk:R.nbins]))
    end
    R.density.solvent = R.solvent.nmols / R.volume.total
    R.density.solute = R.solute.nmols / R.volume.total
    R.density.solvent_bulk = n_solvent_in_bulk / R.volume.bulk

    # Now that we know the the volume of the domain and the density of the solvent in the 
    # bulk region, we can rescale the random counts to take into account that the the ideal
    # gas distribution must have more molecules, with same bulk density, that the true
    # distribution, because we have to take into consieration the available volume which is
    # occupied by the solute
    density_fix = R.density.solvent_bulk / R.density.solvent
    R.md_count_random .= density_fix * R.md_count_random
    R.rdf_count_random .= density_fix * R.rdf_count_random

    # Computing the distribution functions and KB integrals, from the MDDF and from the RDF
    #
    warned_already = false
    for ibin = 1:R.nbins
        if R.md_count_random[ibin] > 0.0
            R.mddf[ibin] = R.md_count[ibin] / R.md_count_random[ibin]
            if ibin == 1
                R.coordination_number[ibin] = R.md_count[ibin]
                R.coordination_number_random[ibin] = R.md_count_random[ibin]
            else
                R.coordination_number[ibin] = R.coordination_number[ibin-1] + R.md_count[ibin]
                R.coordination_number_random[ibin] =
                    R.coordination_number_random[ibin-1] + R.md_count_random[ibin]
            end
        else
            if !warned_already && !options.silent 
                @warn begin
                    """\n
                        Ideal-gas histogram bins with zero samples. 
                        Increase n_random_samples, number of trajectory frames, and/or bin size.

                    """
                end _file=nothing _line=nothing
                warned_already = true
            end
        end
        R.kb[ibin] =
            units.Angs3tocm3permol *
            (1 / R.density.solvent_bulk) *
            (R.coordination_number[ibin] - R.coordination_number_random[ibin])

        # For the RDF
        if R.rdf_count_random[ibin] > 0.0
            R.rdf[ibin] = R.rdf_count[ibin] / R.rdf_count_random[ibin]
        end
        if ibin == 1
            R.sum_rdf_count[ibin] = R.rdf_count[ibin]
            R.sum_rdf_count_random[ibin] = R.rdf_count_random[ibin]
        else
            R.sum_rdf_count[ibin] = R.sum_rdf_count[ibin-1] + R.rdf_count[ibin]
            R.sum_rdf_count_random[ibin] =
                R.sum_rdf_count_random[ibin-1] + R.rdf_count_random[ibin]
        end
        R.kb_rdf[ibin] =
            units.Angs3tocm3permol *
            (1 / R.density.solvent_bulk) *
            (R.sum_rdf_count[ibin] - R.sum_rdf_count_random[ibin])

    end

    return R
end

function _coordination_number_final_results!(R::Result, options::Options)

    if !options.silent
        @warn begin
            """\n
                coordination_number_only was set to true, so the MDDF and KB integrals were not computed. 
                (to remove this warning use `Options(silent=true)`)

            """
        end _file=nothing _line=nothing
    end

    # Setup the distance vector
    R.d .= shellradius.(1:R.nbins, options.binstep)

    # Normalization of number of frames: sum of weights for all frames read
    Q = sum_frame_weights(R)

    # Scale counters by number of samples and frames
    @. R.md_count = R.md_count / (R.solute.nmols * Q)
    @. R.solute_group_count = R.solute_group_count / (R.solute.nmols * Q)
    if R.autocorrelation
        R.solvent_group_count .= R.solute_group_count
    else
        @. R.solvent_group_count = R.solvent_group_count / (R.solute.nmols * Q)
    end
    @. R.rdf_count = R.rdf_count / (R.solute.nmols * Q)

    # Volume of each bin shell and of the solute domain
    R.volume.total = R.volume.total / Q

    R.density.solvent = R.solvent.nmols / R.volume.total
    R.density.solute = R.solute.nmols / R.volume.total
    # Coordination numbers
    for ibin = 1:R.nbins
        if ibin == 1
            R.coordination_number[ibin] = R.md_count[ibin]
            R.sum_rdf_count[ibin] = R.rdf_count[ibin]
        else
            R.coordination_number[ibin] = R.coordination_number[ibin-1] + R.md_count[ibin]
            R.sum_rdf_count[ibin] = R.sum_rdf_count[ibin-1] + R.rdf_count[ibin]
        end
    end

    return R
end

"""
    merge(r::Vector{Result})

This function merges the results of MDDF calculations obtained by running the same
analysis on multiple trajectories, or multiple parts of the same trajectory. It returns
a Result structure of the same type, with all the functions and counters representing averages
of the set provided weighted by the number of frames read in each Result set.

"""
function Base.merge(results::Vector{<:Result})
    cannot_merge = false
    nframes_read = 0
    for ir in eachindex(results)
        for file in results[ir].files
            nframes_read += file.nframes_read
        end
        for jr in ir+1:lastindex(results)
            if results[ir].nbins != results[jr].nbins
                println(
                    "ERROR: To merge Results, the number of bins of the histograms of both sets must be the same.",
                )
                cannot_merge = true
            end
            if !(results[ir].cutoff ≈ results[ir].cutoff)
                println(
                    "ERROR: To merge Results, cutoff distance of the of the histograms of both sets must be the same.",
                )
                cannot_merge = true
            end
            if (results[ir].solute != results[jr].solute) || (results[ir].solvent != results[jr].solvent)
                println(
                    "ERROR: To merge Results, the solute and solvent selections of both sets must be the same.",
                )
                cannot_merge = true
            end
        end
    end
    if cannot_merge
        throw(ArgumentError(" Incompatible set of results to merge. "))
    end

    # Total number of frames of all results, and total weight of frames
    nfiles = 0
    ntot_frames = 0
    tot_frame_weight = 0.0
    for result in results
        nfiles += length(result.files)
        ntot_frames += sum(file.nframes_read for file in result.files)
        tot_frame_weight += sum_frame_weights(result)
    end

    # Compute weight of each file, given the number frames. Also compute 
    # the weight of the data of each file, given the weights of the frames
    weights = zeros(nfiles)
    ifile = 0
    for result in results
        for file in result.files
            ifile += 1
            weights[ifile] = file.nframes_read / ntot_frames
        end
    end

    # Append all TrajectoryFileOptions to single vector
    files = TrajectoryFileOptions[]
    for result in results
        for file in result.files
            push!(files, file)
        end
    end

    # Initialize group counts
    solute_group_count = [ zeros(results[1].nbins) for _ in 1:length(results[1].solute_group_count) ]
    solvent_group_count = [ zeros(results[1].nbins) for _ in 1:length(results[1].solvent_group_count) ] 

    # Structure for merged results
    R = Result(
        nbins = results[1].nbins,
        dbulk = results[1].dbulk,
        cutoff = results[1].cutoff,
        autocorrelation = results[1].autocorrelation,
        solute = results[1].solute,
        solvent = results[1].solvent,
        solute_group_count = solute_group_count,
        solvent_group_count = solvent_group_count,
        files = files,
        weights = weights,
    )

    # Average results weighting the data considering the weights of the frames of each data set
    warn = false
    @. R.d = results[1].d
    for ifile in eachindex(results)
        result = results[ifile]
        w = sum_frame_weights(result) / tot_frame_weight
        if !(w ≈ R.weights[ifile]) && !warn
            warn = true
            @warn begin 
                """\n
                    Frame weights and file weights differ, because crustom frame weights were provided.

                """ 
             end _file=nothing _line=nothing
        end
        @. R.mddf += w * result.mddf
        @. R.kb += w * result.kb
        @. R.rdf += w * result.rdf
        @. R.kb_rdf += w * result.kb_rdf
        @. R.md_count += w * result.md_count
        @. R.md_count_random += w * result.md_count_random
        @. R.coordination_number += w * result.coordination_number
        @. R.coordination_number_random += w * result.coordination_number_random
        for i in eachindex(R.solute_group_count, result.solute_group_count)
            @. R.solute_group_count[i] += w * result.solute_group_count[i]
        end
        for i in eachindex(R.solvent_group_count, result.solvent_group_count)
            @. R.solvent_group_count[i] += w * result.solvent_group_count[i]
        end
        @. R.rdf_count += w * result.rdf_count
        @. R.rdf_count_random += w * result.rdf_count_random
        @. R.sum_rdf_count += w * result.sum_rdf_count
        @. R.sum_rdf_count_random += w * result.sum_rdf_count_random
        R.density.solute += w * result.density.solute
        R.density.solvent += w * result.density.solvent
        R.density.solvent_bulk += w * result.density.solvent_bulk
        R.volume.total += w * result.volume.total
        R.volume.bulk += w * result.volume.bulk
        R.volume.domain += w * result.volume.domain
        R.volume.shell += w * result.volume.shell
    end
    return R
end

@testitem "merge" begin
    using ComplexMixtures: mddf, merge
    using PDBTools: readPDB, select, selindex
    using ComplexMixtures.Testing: data_dir

    # Test simple three-molecule system
    atoms = readPDB("$data_dir/toy/cross.pdb")
    protein = AtomSelection(select(atoms, "protein and model 1"), nmols = 1)
    water = AtomSelection(select(atoms, "resname WAT and model 1"), natomspermol = 3)
    traj = Trajectory("$data_dir/toy/cross.pdb", protein, water, format = "PDBTraj")

    options = Options(
        seed = 321,
        StableRNG = true,
        nthreads = 1,
        silent = true,
        n_random_samples = 10^5,
        lastframe = 1,
    )
    R1 = mddf(traj, options)

    options = Options(
        seed = 321,
        StableRNG = true,
        nthreads = 1,
        silent = true,
        n_random_samples = 10^5,
        firstframe = 2,
    )
    R2 = mddf(traj, options)

    R = merge([R1, R2])
    @test R.volume.total == 27000.0
    @test R.volume.domain ≈ R.volume.total - R.volume.bulk
    @test isapprox(R.volume.domain, (4π / 3) * R.dbulk^3; rtol = 0.01)
    @test R.density.solute ≈ 1 / R.volume.total
    @test R.density.solvent ≈ 3 / R.volume.total
    @test R.density.solvent_bulk ≈ 2 / R.volume.bulk
    @test R.weights == [0.5, 0.5]

    # Test loading a saved merged file
    dir = mktempdir()
    save(R,"$dir/merged.json")
    R_save = load("$dir/merged.json")
    @test isapprox(R, R_save, debug = true)

    # Test merging files for which weights are provided for the frames
    R2 = mddf(traj, options, frame_weights = [0.0, 2.0])
    @test R.weights == [0.5, 0.5]
    @test length(R.files) == 2

    # Two-atom system
    at1 = AtomSelection([1], nmols=1)
    at2 = AtomSelection([2], nmols=1)
    traj = Trajectory("$data_dir/toy/self_monoatomic.pdb", at1, at2, format = "PDBTraj")
    R1 = mddf(traj, Options(lastframe=1))
    @test sum(R1.md_count) == 1
    R2 = mddf(traj, Options(firstframe=2))
    @test sum(R2.md_count) == 0
    R = merge([R1, R2])
    @test R.weights == [0.5, 0.5]
    @test sum(R.md_count) == 0.5
    R1 = mddf(traj, Options(lastframe=1), frame_weights=[2.0])
    @test sum(R1.md_count) == 1
    R = merge([R1, R2])
    @test sum(R.md_count) == 2/3
    @test sum(sum.(R1.solute_group_count)) == 1
    @test sum(sum.(R1.solvent_group_count)) == 1
    @test sum(sum.(R2.solute_group_count)) == 0
    @test sum(sum.(R2.solvent_group_count)) == 0
    @test sum(sum.(R.solute_group_count)) == 2/3
    @test sum(sum.(R.solvent_group_count)) == 2/3

    # Test throwing merging incompatible results
    protein = AtomSelection(select(atoms, "protein and model 1"), nmols = 1)
    traj = Trajectory("$data_dir/toy/cross.pdb", protein, water, format = "PDBTraj")
    R1 = mddf(traj, options)
    protein = AtomSelection(
        select(atoms, "protein and model 1"), nmols = 1, 
        group_names = ["acidic"],
        group_atom_indices = [selindex(atoms, "protein and acidic")]
    )
    traj = Trajectory("$data_dir/toy/cross.pdb", protein, water, format = "PDBTraj")
    R2 = mddf(traj, options)


end

@testitem "Result - empty" begin
    #using ComplexMixtures: Result, Trajectory, Options, AtomSelection
    using ComplexMixtures.Testing: data_dir, pdbfile
    using PDBTools: readPDB, select, name
    atoms = readPDB(pdbfile)
    protein = select(atoms, "protein")
    tmao = select(atoms, "resname TMAO")
    solute = AtomSelection(protein, nmols = 1)
    solvent = AtomSelection(tmao, natomspermol = 14)
    traj = Trajectory("$data_dir/NAMD/trajectory.dcd", solute, solvent)
    options = Options()
    # At this point we can only test an empty Result struct
    R = Result(traj, options)
    @test R.autocorrelation == false
    @test R.cutoff == 10.0
    @test length(R.d) == 500
    @test R.dbulk == 10.0
    @test (R.density.solute, R.density.solvent, R.density.solvent_bulk) == (0.0, 0.0, 0.0)
    @test normpath(R.files[1].filename) == normpath("$data_dir/NAMD/trajectory.dcd")
    @test R.files[1].irefatom == 1
    @test length(R.kb) == 500
    @test length(R.kb_rdf) == 500
    @test R.files[1].lastframe_read == 20
    @test length(R.md_count) == 500
    @test length(R.md_count_random) == 500
    @test length(R.mddf) == 500
    @test R.nbins == 500
    @test R.files[1].nframes_read == 20
    @test R.files[1].options == Options()
    @test length(R.rdf_count) == 500
    @test length(R.rdf_count_random) == 500
    @test R.solute == AtomSelection(collect(1:1463), nmols=1, natomspermol=1463, group_names=name.(protein))
    @test R.solvent == AtomSelection(collect(1479:4012), nmols=181, natomspermol=14, group_names=name.(tmao[1:14]))
    @test length(R.solute_group_count) == 1463
    @test all(length.(R.solute_group_count) .== 500)
    @test length(R.solvent_group_count) == 14
    @test all(length.(R.solvent_group_count) .== 500)
    @test length(R.coordination_number) == 500
    @test length(R.coordination_number_random) == 500
    @test length(R.sum_rdf_count_random) == 500
    @test (R.volume.total, R.volume.bulk, R.volume.domain, length(R.volume.shell)) == (0.0, 0.0, 0.0, 500)
    @test R.weights[1] == 1.0
end

#
# Functions to save the results to a file
#
# These definitions are probably not necessary for StructTypes >= 1.10.0
# 
StructTypes.StructType(::Type{AtomSelection}) = StructTypes.Struct()
StructTypes.StructType(::Type{Result}) = StructTypes.Struct()
StructTypes.StructType(::Type{Density}) = StructTypes.Struct()
StructTypes.StructType(::Type{Volume}) = StructTypes.Struct()
StructTypes.StructType(::Type{Options}) = StructTypes.Struct()
StructTypes.StructType(::Type{TrajectoryFileOptions}) = StructTypes.Struct()

"""
    save(R::Result, filename::String)

Function to write the result data structure to a json file.

"""
function save(filename::String, R::Result)
    filename = expanduser(filename)
    open(filename, "w") do f
        JSON3.write(f, R)
    end
    return "Results saved in JSON file: $filename"
end
# legacy order
save(R::Result, filename::String) = save(filename, R)

#
# This function tries to read a version number from a result.json
# file. This is used to adjust the reading to legacy formats of the
# output file.
#
function _get_version(filename)
    str = readuntil(expanduser(filename), ',')
    v = match(r"\"Version\":\"([^\"]*)\"", str)
    return isnothing(v) ? v"1.0.0" : VersionNumber(v[1])
end

"""
    load(filename::String, ::Type{T}=Result)

Function to load the json saved results file into, by default, the `Result` data structure.

## Example

```julia
using ComplexMixtures: load
R = load("results.json")
#or
R = load("results.json", Result)
```

"""
function load(filename::String, ::Type{Result}=Result)
    filename = expanduser(filename)
    json_version = _get_version(filename)
    current_version = pkgversion(@__MODULE__)
    # Error if the json file is from a newer version than the current one
    if json_version > current_version
        throw(ArgumentError("""\n
            Trying to load a json result file created with a newer version of ComplexMixtures. 
            This can cause unpredictable errors. 

            Current version of ComplexMixtures: $current_version
            Version used to create the output .json file: $json_version

            Please update ComplexMixtures and try again.

        """))
    end
    if json_version < v"2.0.0"
        throw(ArgumentError("""\n
            Trying to load a json result created with an older, and incompatible, version of ComplexMixtures.

            Current version of ComplexMixtures: $current_version
            Version used to create the output .json file: $json_version

            To load this file, install an older version of ComplexMixtures, with, for example:
            
            julia> import Pkg; Pkg.pkg"add ComplexMixtures@$json_version"

            You can pin the version of ComplexMixtures to the one you installed with:

            julia> import Pkg; Pkg.pkg"pin ComplexMixtures@$json_version"
        
        """))
    end
    R = open(filename, "r") do io
        JSON3.read(io, Result)
    end
    return R
end

@testitem "Result - load/save" begin
    using ComplexMixtures: load
    using ComplexMixtures.Testing: data_dir
    r1 = load("$data_dir/NAMD/protein_tmao.json")
    tmp = tempname()
    save(r1, tmp)
    r2 = load(tmp)
    @test r1 == r2
    r2 = load(tmp, Result)
    @test r1 == r2
    save(tmp, r1)
    r2 = load(tmp)
    @test r1 == r2
    # Test throwing an error incompatible versions of ComplexMixtures
    @test_throws ArgumentError load("$data_dir/wrong_version_jsons/too_new.json")
    @test_throws ArgumentError load("$data_dir/wrong_version_jsons/too_old.json")
end


#=
    which_types(s::AtomSelection, indices::Vector{Int})

Function that returns the list of the indices of the types of the atoms in a
selection. For example, if a selection corresponds to a solvent of water molecules:
There are three types, 1, 2, and 3, corresponding to the three atoms of the
water molecule. If the indices provided are, for instance, 11, 12, and 13, 
corresponding to a water molecule, this function will return 1, 2 and 3.

This is used to get equivalent-atom contributions to the distribution functions.
For example, the input indices span all water molecules, the output of this
function will be still the three indices corresponding to the three types
of atoms that exist in a water molecule. 

It is not possible to compute the contribution of *one* individual water molecule
if the distribution function was computed for all molecules. Thus, the necessity
to identify the types of atoms involved in a selection.   

=#
function which_types(s::AtomSelection, indices::Vector{Int}; warning = true)
    selected_types = Int[]
    ntypes = 0
    for i in indices
        isel = findfirst(ind -> ind == i, s.indices)
        if isnothing(isel)
            error(" Atom in input list is not part of solvent (or solute).")
        else
            it = atom_type(isel, s.natomspermol)
            if !(it in selected_types)
                push!(selected_types, it)
                ntypes += 1
                if ntypes == s.natomspermol
                    warning && println(
                        "WARNING: All possible types of atoms ($ntypes) of this selection were selected.",
                    )
                    return selected_types
                end
            end
        end
    end
    return selected_types
end

#=
    sum!(R1::Result, R2::Result)

Sum the counts of two Results structures, adding the result to the first structure as in R1 = R1 + R2.

=#
function sum!(R1::Result, R2::Result)

    @. R1.md_count += R2.md_count
    @. R1.md_count_random += R2.md_count_random

    for i in eachindex(R1.solute_group_count)
        @. R1.solute_group_count[i] += R2.solute_group_count[i]
    end
    for i in eachindex(R1.solvent_group_count)
        @. R1.solvent_group_count[i] += R2.solvent_group_count[i]
    end

    @. R1.rdf_count += R2.rdf_count
    @. R1.rdf_count_random += R2.rdf_count_random

    R1.volume.total += R2.volume.total

    return R1
end

#=
    title(R::Result, solute::AtomSelection, solvent::AtomSelection)
    title(R::Result, solute::AtomSelection, solvent::AtomSelection, nspawn::Int)

Print some information about the run.

=#
function title(R::Result, solute::AtomSelection, solvent::AtomSelection)
    print(
        """
        $(bars)
        Starting MDDF calculation:
        $(R.files[1].nframes_read) frames will be considered.
        Solute: $(atoms_str(natoms(solute))) belonging to $(mol_str(solute.nmols)).
        Solvent: $(atoms_str(natoms(natoms))) belonging to $(mol_str(solvent.nmols))
        """)
end
function title(R::Result, solute::AtomSelection, solvent::AtomSelection, nspawn::Int)
    print(
        """ 
        $(bars)
        Starting MDDF calculation in parallel:
        $(R.files[1].nframes_read) frames will be considered.
        Number of calculation threads: $(nspawn)
        Solute: $(atoms_str(natoms(solute))) belonging to $(mol_str(solute.nmols)).
        Solvent: $(atoms_str(natoms(solvent))) belonging to $(mol_str(solvent.nmols)).
        """)
end

#
# Print overview of the results in the REPL
#
#=

Structure that is used to dispatch the show of a overview.

=#
@kwdef mutable struct Overview
    R::Result
    domain_molar_volume::Float64 = 0.0
    density::Density = Density()
    solvent_molar_volume::Float64 = 0.0
    solvent_molar_volume_bulk::Float64 = 0.0
    solute_molar_volume::Float64 = 0.0
end

function _bulk_range_from_R(R)
    if R.dbulk == R.cutoff
        return ">= $(R.dbulk) Å"
    else
        return "$(R.dbulk) - $(R.cutoff) Å"
    end
end

function Base.show(io::IO, ov::Overview)
    println(
        io,
        """
 $bars
 MDDF Overview - ComplexMixtures - Version $(ov.R.Version)
 $bars

 Solvent properties:
 -------------------

 Simulation concentration: $(ov.density.solvent) mol L⁻¹
 Molar volume: $(ov.solvent_molar_volume) cm³ mol⁻¹

 Concentration in bulk: $(ov.density.solvent_bulk) mol L⁻¹
 Molar volume in bulk: $(ov.solvent_molar_volume_bulk) cm³ mol⁻¹

 Solute properties:
 ------------------

 Simulation Concentration: $(ov.density.solute) mol L⁻¹
 Estimated solute partial molar volume: $(ov.solute_molar_volume) cm³ mol⁻¹

 Bulk range: $(_bulk_range_from_R(ov.R))
 Molar volume of the solute domain: $(ov.domain_molar_volume) cm³ mol⁻¹

 Auto-correlation: $(ov.R.autocorrelation)

 Trajectory files and weights:
 """
    )
    for i in eachindex(ov.R.files)
        println(io, "   $(normpath(ov.R.files[i].filename)) - w = $(ov.R.weights[i])")
    end
    ifar = trunc(Int, ov.R.nbins - 1.0 / ov.R.files[1].options.binstep)
    long_range_mean = mean(ov.R.mddf[ifar:ov.R.nbins])
    long_range_std = std(ov.R.mddf[ifar:ov.R.nbins])
    long_range_mean_rdf = mean(ov.R.rdf[ifar:ov.R.nbins])
    long_range_std_rdf = std(ov.R.rdf[ifar:ov.R.nbins])
    print(
        io,
        """

   Long range MDDF mean (expected 1.0): $long_range_mean ± $long_range_std
   Long range RDF mean (expected 1.0): $long_range_mean_rdf ± $long_range_std_rdf

   $bars"""
   )
end

"""
    overview(R::Result)

Function that outputs the volumes and densities in the most natural units.
"""
function overview(R::Result)

    ov = Overview(R = R)

    # Molar volume of the solute domain
    ov.domain_molar_volume = R.volume.domain * units.Angs3tocm3permol

    # Density of the solute and of the solvent 
    ov.density.solute = R.density.solute * units.SitesperAngs3tomolperL
    ov.density.solvent = R.density.solvent * units.SitesperAngs3tomolperL
    ov.density.solvent_bulk = R.density.solvent_bulk * units.SitesperAngs3tomolperL

    # Solvent molar volume ov.solvent_molar_volume = 1000 / ov.density.solvent
    ov.solvent_molar_volume = 1000 / ov.density.solvent
    ov.solvent_molar_volume_bulk = 1000 / ov.density.solvent_bulk

    # Solute molar volume computed from solvent density in bulk
    if R.autocorrelation
        ov.solute_molar_volume = ov.solvent_molar_volume
    else
        ov.solute_molar_volume =
            units.Angs3tocm3permol *
            (R.density.solvent_bulk * R.volume.total - R.solvent.nmols) /
            R.density.solvent_bulk
    end

    return ov
end
