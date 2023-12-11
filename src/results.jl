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

"""
    setbin(d,step)

$(INTERNAL)

Function that sets to which histogram bin a data point pertains simple, but important to keep consistency over all calls.

"""
setbin(d, step) = max(1, ceil(Int, d / step))

"""

$(TYPEDEF)

Structures to contain the volumes obtained from calculations.

$(TYPEDFIELDS)

"""
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

$(INTERNAL)

$(TYPEDEF)

Structures to contain the details of a solute or solvent to store in the results of the MDDF calculation.

$(TYPEDFIELDS)

"""
struct SolSummary
    natoms::Int
    nmols::Int
    natomspermol::Int
end
SolSummary(s::Selection) = SolSummary(s.natoms, s.nmols, s.natomspermol)

function Base.show(io::IO, s::SolSummary)
    println(io, "Number of atoms: ", s.natoms)
    println(io, "Number of molecules: ", s.nmols)
    print(io, "Number of atoms per molecule: ", s.natomspermol)
end

"""

$(TYPEDEF)

Structure to contain the results of the MDDF calculation.

$(TYPEDFIELDS)

The Result{Vector{Float64}} parametric type is necessary only for reading the JSON3 saved file. 

"""
@kwdef mutable struct Result{T<:VecOrMat{Float64}}

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
    mddf::Vector{Float64} = zeros(nbins)
    kb::Vector{Float64} = zeros(nbins)

    # Properties of the solute and solvent selections
    autocorrelation::Bool
    solvent::SolSummary
    solute::SolSummary

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
    density::Density = Density()
    volume::Volume = Volume(nbins)

    # Options of the calculation
    options::Options
    irefatom::Int
    lastframe_read::Int
    nframes_read::Int

    # File name(s) of the trajectories in this results 
    files::Vector{String}
    weights::Vector{Float64}
end

#
# Initialize the data structure that is returned from the computation, and checks some
# input parameters for consistency
#
function Result(trajectory::Trajectory, options::Options; irefatom = -1)

    # Check for simple input errors
    if options.stride < 1
        error("in MDDF options: stride cannot be less than 1. ")
    end
    if options.lastframe > 0 && options.lastframe < options.firstframe
        error("in MDDF options: lastframe must be greater or equal to firstframe. ")
    end
    if options.lastframe > trajectory.nframes
        error("in MDDF options: lastframe is greater than trajectory.nframes. ")
    end

    # Check for problems in dbulk and cutoff definitions
    cutoff = options.dbulk
    if (options.dbulk / options.binstep) % 1 > 1.e-5
        error("in MDDF options: dbulk must be a multiple of binstep.")
    end
    if options.usecutoff
        if options.dbulk >= options.cutoff
            error(" in MDDF options: The bulk volume is zero (dbulk must be smaller than cutoff). ")
        end
        if (options.cutoff / options.binstep) % 1 > 1.e-5
            error("in MDDF options: cutoff must be a multiple of binstep.")
        end
        if ((options.cutoff - options.dbulk) / options.binstep) % 1 > 1.e-5
            error("in MDDF options: (cutoff-dbulk) must be a multiple of binstep. ")
        end
        cutoff = options.cutoff
    end
    nbins = setbin(cutoff, options.binstep)

    if options.irefatom > trajectory.solvent.natoms
        error("in MDDF options: Reference atom index", options.irefatom, " is greater than number of atoms of the solvent molecule. ")
    end

    # Open trajectory to read some data
    opentraj!(trajectory)
    firstframe!(trajectory)

    # Set reference atom as the closest one to the center of coordinates of the molecule, as default
    if options.irefatom == -1
        nextframe!(trajectory)
        first_mol = viewmol(1, trajectory.x_solvent, trajectory.solvent)
        # don't use findmin for compatibility with Julia 1.6
        irefatom = 0
        dmin = +Inf
        cm = mean(first_mol)
        for (i, at) in pairs(first_mol)
            d = norm(at - cm)
            if d < dmin 
                dmin = d
                irefatom = i
            end
        end
    else
        irefatom = options.irefatom
    end

    # Last frame to be considered
    if options.lastframe == -1
        lastframe_read = trajectory.nframes
    else
        lastframe_read = options.lastframe
    end

    # Close trajecotory
    closetraj!(trajectory)

    # Actual number of frames that are read considering lastframe and stride
    nframes_read = round(Int, (lastframe_read - options.firstframe + 1) / options.stride)
    if nframes_read == 0
        error("Number of frames to read is zero. Check input parameters.")
    end

    return Result(
        options = options,
        nbins = nbins,
        dbulk = options.dbulk,
        cutoff = cutoff,
        irefatom = irefatom,
        lastframe_read = lastframe_read,
        nframes_read = nframes_read,
        autocorrelation = isautocorrelation(trajectory),
        solute = SolSummary(trajectory.solute),
        solvent = SolSummary(trajectory.solvent),
        files = [trajectory.filename],
        weights = [1.0],
    )
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
        samples = (solvent_nmols = R.solvent.nmols - 1, random = R.options.n_random_samples)
    else
        samples = (solvent_nmols = R.solvent.nmols, random = R.options.n_random_samples)
    end
    return samples
end

# autocorrelation can be obtained from the comparison of solute and solvent indexes
isautocorrelation(solute_indexes::Vector{Int}, solvent_indexes::Vector{Int}) =
    solute_indexes == solvent_indexes ? true : false
isautocorrelation(trajectory::Trajectory) =
    isautocorrelation(trajectory.solute.index, trajectory.solvent.index)

#
# Functions to compute volumes of shells
#
"""
    sphericalshellvolume(i,step)

$(INTERNAL)

Computes the volume of the spherical shell defined within [(i-1)*step,i*step].

"""
function sphericalshellvolume(i, step)
    rmin = (i - 1) * step
    return (4 * pi / 3) * ((rmin + step)^3 - rmin^3)
end

"""
    shellradius(i,step)

$(INTERNAL)

Compute the point in which the radius comprises half of the volume of the shell.

"""
function shellradius(i, step)
    rmin = (i - 1) * step
    return (0.5 * ((rmin + step)^3 + rmin^3))^(1 / 3)
end


"""
    sphereradiusfromshellvolume(volume,step)

$(INTERNAL)

Computes the radius that corresponds to a spherical shell of a given volume.

"""
function sphereradiusfromshellvolume(volume, step)
    fourthirdsofpi = 4 * pi / 3
    if 3 * step * volume - pi * step^4 <= 0.0
        return 0.0
    end
    rmin = (sqrt(3 * pi) * sqrt(3 * step * volume - pi * step^4) - 3 * pi * step^2) / (6 * pi * step)
    return (0.5 * (volume / fourthirdsofpi + 2 * rmin^3))^(1 / 3)
end

#
# Function to compute the normalization of the weights, given the optional
# frame weights and the number of frames read, which is used in the case
# that frame weights were not provided
#
function sum_frame_weights(R::Result)
    i = R.options.firstframe
    s = R.options.stride
    l = R.lastframe_read
    Q = if !isempty(R.options.frame_weights)
        sum(R.options.frame_weights[i] for i in i:s:l)
    else
        # number of frames that were red from the file
        round(Int, (l - i + 1) / s)
    end
    return Q
end

"""
    finalresults!(R::Result, options::Options, trajectory::Trajectory)

$(INTERNAL)

Function that computes the final results of all the data computed by averaging according to the sampling of each type of data, and converts to common units.

Computes also the final distribution functions and KB integrals.

This function modified the values contained in the R data structure.

"""
function finalresults!(R::Result, options::Options, trajectory::Trajectory)

    # Sampling scheme depending on the type of calculation
    samples = set_samples(R)

    # Setup the distance vector
    for i = 1:R.nbins
        R.d[i] = shellradius(i, options.binstep)
    end

    # Normalization of number of frames: sum of weights for all frames read
    Q = sum_frame_weights(R)

    # Scale counters by number of samples and frames
    @. R.md_count = R.md_count / (R.solute.nmols * Q)
    @. R.solute_atom = R.solute_atom / (R.solute.nmols * Q)
    @. R.solvent_atom = R.solvent_atom / (R.solute.nmols * Q)
    @. R.md_count_random = R.md_count_random / (samples.random * Q)
    @. R.rdf_count = R.rdf_count / (R.solute.nmols * Q)
    @. R.rdf_count_random = R.rdf_count_random / (samples.random * Q)

    # Volume of each bin shell and of the solute domain
    R.volume.total = R.volume.total / Q
    @. R.volume.shell = R.volume.total * (R.rdf_count_random / samples.solvent_nmols)

    # Solute domain volume
    ibulk = setbin(R.dbulk + 0.5 * R.options.binstep, R.options.binstep)
    R.volume.domain = sum(@view(R.volume.shell[1:ibulk-1]))

    # Bulk volume and density properties: either the bulk is considered everything
    # that is not the domain, or the bulk is the region between d_bulk and cutoff,
    # if R.options.usecutoff is true (meaning that there is a cutoff different from
    # that of the bulk distance)
    if !R.options.usecutoff
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

    #
    # Computing the distribution functions and KB integrals, from the MDDF and from the RDF
    #
    warn = false
    for ibin = 1:R.nbins
        # For the MDDF
        if R.md_count_random[ibin] > 0.0
            R.mddf[ibin] = R.md_count[ibin] / R.md_count_random[ibin]
            for i = 1:trajectory.solute.natomspermol
                R.solute_atom[ibin, i] = R.solute_atom[ibin, i] / R.md_count_random[ibin]
            end
            for j = 1:trajectory.solvent.natomspermol
                R.solvent_atom[ibin, j] = R.solvent_atom[ibin, j] / R.md_count_random[ibin]
            end
        else
            if !warn && !options.silent
                @warn begin
                    """
                    Ideal-gas histogram bins with zero samples. 
                    Increase n_random_samples or trajectory length.
                    """
                end _file=nothing _line=nothing
                warn = true
            end
        end
        if ibin == 1
            R.coordination_number[ibin] = R.md_count[ibin]
            R.coordination_number_random[ibin] = R.md_count_random[ibin]
        else
            R.coordination_number[ibin] = R.coordination_number[ibin-1] + R.md_count[ibin]
            R.coordination_number_random[ibin] =
                R.coordination_number_random[ibin-1] + R.md_count_random[ibin]
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

"""
    merge(r::Vector{Result})

This function merges the results of MDDF calculations obtained by running the same
analysis on multiple trajectories, or multiple parts of the same trajectory. It returns
a Result structure of the same type, with all the functions and counters representing averages
of the set provided weighted by the number of frames read in each Result set.

"""
function Base.merge(r::Vector{<:Result})
    nr = length(r)
    nframes_read = r[1].nframes_read
    error = false
    for ir = 2:nr
        nframes_read += r[ir].nframes_read
        if r[ir].nbins != r[1].nbins
            println(
                "ERROR: To merge Results, the number of bins of the histograms of both sets must be the same.",
            )
        end
        if !(r[ir].cutoff ≈ r[1].cutoff)
            println(
                "ERROR: To merge Results, cutoff distance of the of the histograms of both sets must be the same.",
            )
        end
    end
    if error
        error(" Incompatible set of results to merge. ")
    end

    # List of files and weights
    nfiles = sum(length(R.files) for R in r)
    files = Vector{String}(undef, nfiles)
    weights = Vector{Float64}(undef, nfiles)

    # First, merge the options
    options = merge(getfield.(r, :options))

    # Structure for merged results
    R = Result(
        options = options,
        nbins = r[1].nbins,
        dbulk = r[1].dbulk,
        cutoff = r[1].cutoff,
        irefatom = r[1].irefatom,
        lastframe_read = r[nr].lastframe_read,
        nframes_read = nframes_read,
        autocorrelation = r[1].autocorrelation,
        solute = r[1].solute,
        solvent = r[1].solvent,
        files = files,
        weights = weights,
    )

    # Total normalization factor: sum of the number of frame reads,
    # or the sum of frame weights
    Q = sum_frame_weights(R)

    # Average results weighting the data considering the weights of the frames of each data set
    @. R.d = r[1].d
    ifile = 0
    for ir = 1:nr
        w = sum_frame_weights(r[ir]) / Q
        @. R.mddf += w * r[ir].mddf
        @. R.kb += w * r[ir].kb
        @. R.rdf += w * r[ir].rdf
        @. R.kb_rdf += w * r[ir].kb_rdf
        @. R.md_count += w * r[ir].md_count
        @. R.md_count_random += w * r[ir].md_count_random
        @. R.coordination_number += w * r[ir].coordination_number
        @. R.coordination_number_random += w * r[ir].coordination_number_random
        @. R.solute_atom += w * r[ir].solute_atom
        @. R.solvent_atom += w * r[ir].solvent_atom
        @. R.rdf_count += w * r[ir].rdf_count
        @. R.rdf_count_random += w * r[ir].rdf_count_random
        @. R.sum_rdf_count += w * r[ir].sum_rdf_count
        @. R.sum_rdf_count_random += w * r[ir].sum_rdf_count_random
        R.density.solute += w * r[ir].density.solute
        R.density.solvent += w * r[ir].density.solvent
        R.density.solvent_bulk += w * r[ir].density.solvent_bulk
        R.volume.total += w * r[ir].volume.total
        R.volume.bulk += w * r[ir].volume.bulk
        R.volume.domain += w * r[ir].volume.domain
        R.volume.shell += w * r[ir].volume.shell
        for j = 1:length(r[ir].files)
            ifile += 1
            R.files[ifile] = normpath(r[ir].files[j])
            R.weights[ifile] = w * r[ir].weights[j]
        end
    end
    return R
end

@testitem "merge" begin
    using ComplexMixtures
    using PDBTools
    using ComplexMixtures.Testing

    # Test simple three-molecule system
    atoms = readPDB("$(Testing.data_dir)/toy/cross.pdb")
    protein = Selection(select(atoms, "protein and model 1"), nmols = 1)
    water = Selection(select(atoms, "resname WAT and model 1"), natomspermol = 3)
    traj = Trajectory("$(Testing.data_dir)/toy/cross.pdb", protein, water, format = "PDBTraj")

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

    # Test loading a saved merged file
    dir = "$(Testing.data_dir)/NAMD"
    atoms = readPDB("$dir/structure.pdb")
    tmao = Selection(select(atoms, "resname TMAO"), natomspermol = 14)
    water = Selection(select(atoms, "water"), natomspermol = 3)

    #save(R,"$dir/merged.json")
    R_save = load("$dir/merged.json")

    options = Options(
        firstframe = 1,
        lastframe = 2,
        seed = 321,
        StableRNG = true,
        nthreads = 1,
        silent = true,
    )
    traj = Trajectory("$dir/trajectory.dcd", tmao, water)
    R1 = mddf(traj, options)

    options = Options(
        firstframe = 3,
        lastframe = 6,
        seed = 321,
        StableRNG = true,
        nthreads = 1,
        silent = true,
    )
    traj = Trajectory("$dir/trajectory.dcd", tmao, water)
    R2 = mddf(traj, options)

    R = merge([R1, R2])
    @test isapprox(R, R_save, debug = true)

    # Test merging files for which weights are provided for the frames
    traj = Trajectory("$dir/trajectory.dcd", tmao, water)
    options = Options(firstframe = 1, lastframe = 2, seed = 321, StableRNG = true, nthreads = 1, silent = true, frame_weights = fill(0.3, 20))
    R1 = mddf(traj, options)

    # First lets test the error message, in case the frame_weights are not provided for all frames
    options = Options(firstframe = 1, lastframe = 2, seed = 321, StableRNG = true, nthreads = 1, silent = true)
    R2 = mddf(traj, options)
    @test_throws ArgumentError merge([R1, R2])

end

@testitem "Result - empty" begin
    using ComplexMixtures
    using ComplexMixtures.Testing
    using PDBTools
    atoms = readPDB(Testing.pdbfile)
    protein = Selection(select(atoms, "protein"), nmols = 1)
    tmao = Selection(select(atoms, "resname TMAO"), natomspermol = 14)
    traj = Trajectory("$(Testing.data_dir)/NAMD/trajectory.dcd", protein, tmao)
    options = Options()
    # At this point we can only test an empty Result struct
    R = Result(traj, options)
    @test R.autocorrelation == false
    @test R.cutoff == 10.0
    @test length(R.d) == 500
    @test R.dbulk == 10.0
    @test (R.density.solute, R.density.solvent, R.density.solvent_bulk) == (0.0, 0.0, 0.0)
    @test normpath(R.files[1]) == normpath("$(Testing.data_dir)/NAMD/trajectory.dcd")
    @test R.irefatom == 1
    @test length(R.kb) == 500
    @test length(R.kb_rdf) == 500
    @test R.lastframe_read == 20
    @test length(R.md_count) == 500
    @test length(R.md_count_random) == 500
    @test length(R.mddf) == 500
    @test R.nbins == 500
    @test R.nframes_read == 20
    @test R.options == Options()
    @test length(R.rdf_count) == 500
    @test length(R.rdf_count_random) == 500
    @test R.solute == ComplexMixtures.SolSummary(1463, 1, 1463)
    @test R.solvent == ComplexMixtures.SolSummary(2534, 181, 14)
    @test size(R.solute_atom) == (500, 1463)
    @test size(R.solvent_atom) == (500, 14)
    @test length(R.coordination_number) == 500
    @test length(R.coordination_number_random) == 500
    @test length(R.sum_rdf_count_random) == 500
    @test (R.volume.total, R.volume.bulk, R.volume.domain, length(R.volume.shell)) == (0.0, 0.0, 0.0, 500)
    @test R.weights[1] == 1.0
end

#
# Functions to save the results to a file
#
StructTypes.StructType(::Type{SolSummary}) = StructTypes.Struct()
StructTypes.StructType(::Type{Result{Vector{Float64}}}) = StructTypes.Struct()
StructTypes.StructType(::Type{Result{Matrix{Float64}}}) = StructTypes.Struct()
StructTypes.StructType(::Type{Density}) = StructTypes.Struct()
StructTypes.StructType(::Type{Volume}) = StructTypes.Struct()
StructTypes.StructType(::Type{Options}) = StructTypes.Struct()

"""
    save(R::Result, filename::String)

Function to write the result data structure to a json file.

"""
function save(R::Result, filename::String)
    open(filename, "w") do f
        JSON3.write(f, R)
    end
    return "Results saved in JSON file: $filename"
end

#
# This function tries to read a version number from a result.json
# file. This is used to adjust the reading to legacy formats of the
# output file.
#
function _get_version(filename)
    str = readuntil(filename, ',')
    v = match(r"\"Version\":\"([^\"]*)\"", str)
    return isnothing(v) ? v"1.0.0" : VersionNumber(v[1])
end

"""
    load(filename::String)

Function to load the json saved results file into the `Result` data structure.

"""
function load(filename::String; legacy_warning = true)
    json_version = _get_version(filename)
    current_version = pkgversion(@__MODULE__)
    # Error if the json file is from a newer version than the current one
    if json_version > current_version
        error("""
            Trying to load a json result file created with a newer version of ComplexMixtures. 
            This can cause unpredictable errors. 

            Current version of ComplexMixtures: $current_version
            Version used to create the output .json file: $json_version

            Please update ComplexMixtures and try again.
        """)
    end
    # Load directly if within the output compatibility threshold 
    if json_version > v"1.3.4"
        f = open(filename, "r")
        R = JSON3.read(f, Result{Vector{Float64}})
        close(f)
        # Need to reshape the solute and solvent atom contributions, because the data is read in a single column
        solute_atom = reshape(R.solute_atom, R.nbins, :)
        solvent_atom = reshape(R.solvent_atom, R.nbins, :)
        r_names = fieldnames(Result)
        # Return the Result{Matrix{Float64}} type with the appropriate fields 
        return Result{Matrix{Float64}}(
            ntuple(length(r_names)) do i
                r_names[i] == :solute_atom ? solute_atom :
                r_names[i] == :solvent_atom ? solvent_atom : getfield(R, r_names[i])
            end...,
        )
    end
    # Load legacy results
    json_version_str = json_version >= v"1.3.5" ? json_version : "<= 1.3.4"
    if legacy_warning
        @warn begin 
            """\n
            LOADING RESULT JSON FILE IN LEGACY FORMAT. 
    
            Current version of ComplexMixtures: $current_version
            Version used to create the json file: $json_version_str
    
            If the current version overwrites the json file, it may not be readable
            with the older version of ComplexMixtures used to originally create it. 
    
            Note: Output files generated with versions older than 1.0.0 will error.
    
            You can disable this warning by using `load(filename; legacy_warning = false)`
    
            """
        end _file=nothing _line=nothing
    end
    results_updated = load_legacy_json(filename, json_version) 
    return results_updated
end

@testitem "Result - load/save" begin
    using ComplexMixtures
    using ComplexMixtures.Testing
    r1 = load("$(Testing.data_dir)/NAMD/protein_tmao.json", legacy_warning = false)
    tmp = tempname()
    save(r1, tmp)
    r2 = load(tmp)
    @test r1 == r2
end


"""
    which_types(s::Selection, indexes::Vector{Int})

$(INTERNAL)

Function that returns the list of the indexes of the types of the atoms in a
selection. For example, if a selection corresponds to a solvent of water molecules:
There are three types, 1, 2, and 3, corresponding to the three atoms of the
water molecule. If the indexes provided are, for instance, 11, 12, and 13, 
corresponding to a water molecule, this function will return 1, 2 and 3.

This is used to get equivalent-atom contributions to the distribution functions.
For example, the input indexes span all water molecules, the output of this
function will be still the three indexes corresponding to the three types
of atoms that exist in a water molecule. 

It is not possible to compute the contribution of *one* individual water molecule
if the distribution function was computed for all molecules. Thus, the necessity
to identify the types of atoms involved in a selection.   

"""
function which_types(s::Selection, indexes::Vector{Int}; warning = true)
    selected_types = Vector{Int}(undef, 0)
    ntypes = 0
    for i in indexes
        isel = findfirst(ind -> ind == i, s.index)
        if isnothing(isel)
            error(" Atom in input list is not part of solvent (or solute).")
        else
            it = itype(isel, s.natomspermol)
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

"""
    sum!(R1::Result, R2::Result)

$(INTERNAL)

Sum the counts of two Results structures, adding the result to the first structure as in R1 = R1 + R2.

"""
function sum!(R1::Result, R2::Result)

    @. R1.md_count += R2.md_count
    @. R1.md_count_random += R2.md_count_random

    @. R1.solute_atom += R2.solute_atom
    @. R1.solvent_atom += R2.solvent_atom

    @. R1.rdf_count += R2.rdf_count
    @. R1.rdf_count_random += R2.rdf_count_random

    R1.volume.total += R2.volume.total

    return R1
end

"""
    title(R::Result, solute::Selection, solvent::Selection)
    title(R::Result, solute::Selection, solvent::Selection, nspawn::Int)

$(INTERNAL)

Print some information about the run.

"""
function title(R::Result, solute::Selection, solvent::Selection)
    print(
        """
        $(bars)
        Starting MDDF calculation:
        $(R.nframes_read) frames will be considered.
        Solute: $(atoms_str(solute.natoms)) belonging to $(mol_str(solute.nmols)).
        Solvent: $(atoms_str(solvent.natoms)) belonging to $(mol_str(solvent.nmols))
        """)
end
function title(R::Result, solute::Selection, solvent::Selection, nspawn::Int)
    print(
        """ 
        $(bars)
        Starting MDDF calculation in parallel:
        $(R.nframes_read) frames will be considered.
        Number of calculation threads: $(nspawn)
        Solute: $(atoms_str(solute.natoms)) belonging to $(mol_str(solute.nmols)).
        Solvent: $(atoms_str(solvent.natoms)) belonging to $(mol_str(solvent.nmols)).
        """)
end

#
# Print overview of the results in the REPL
#
"""

$(INTERNAL)

$(TYPEDEF)

Structure that is used to dispatch the show of a overview.

$(TYPEDFIELDS)

"""
@kwdef mutable struct Overview
    R::Result
    domain_molar_volume::Float64 = 0.0
    density::Density = Density()
    solvent_molar_volume::Float64 = 0.0
    solvent_molar_volume_bulk::Float64 = 0.0
    solute_molar_volume::Float64 = 0.0
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

 Using dbulk = $(ov.R.dbulk)Å:
 Molar volume of the solute domain: $(ov.domain_molar_volume) cm³ mol⁻¹

 Auto-correlation: $(ov.R.autocorrelation)

 Trajectory files and weights:
 """
    )
    for i = 1:length(ov.R.files)
        println(io, "   $(normpath(ov.R.files[i])) - w = $(ov.R.weights[i])")
    end
    ifar = trunc(Int, ov.R.nbins - 1.0 / ov.R.options.binstep)
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

    # Solvent molar volume
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
