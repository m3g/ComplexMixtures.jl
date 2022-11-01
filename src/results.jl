"""

$(TYPEDEF)

Structure to contain the density values obtained from the calculation.

$(TYPEDFIELDS)

"""
@with_kw mutable struct Density
    solute::Float64 = 0.0
    solvent::Float64 = 0.0
    solvent_bulk::Float64 = 0.0
end

function reset!(d::Density)
    d.solute = 0.0
    d.solvent = 0.0
    d.solvent_bulk = 0.0
    return nothing
end

#function Base.show(io::IO, d::Density ) 
#  println(" Mean solute density: $(d.solute) ")
#  println(" Mean solvent density: $(d.solvent) ")
#  println(" Mean solvent bulk density: $(d.solvent_bulk) ")
#end


"""

$(TYPEDEF)

Structures to contain the volumes obtained from calculations.

$(TYPEDFIELDS)

"""
@with_kw mutable struct Volume
    total::Float64
    bulk::Float64
    domain::Float64
    shell::Vector{Float64}
end

Volume(nbins::Int) = Volume(0.0, 0.0, 0.0, zeros(Float64, nbins))

function reset!(v::Volume)
    v.total = 0.0
    v.bulk = 0.0
    v.domain = 0.0
    @. v.shell = 0.0
    return nothing
end

#function Base.show(io::IO, v::Volume) 
#  n = length(v.shell)
#  println(" Mean total box volume: $(v.total) ")
#  println(" Mean bulk volume: $(v.bulk) ")
#  println(" Mean solute domain volume: $(v.domain) ")
#  println(" Volumes of first, medium, and last solvation shells: $(v.shell[1]), $(v.shell[round(Int,n/2)]), $(v.shell[n])")
#end


"""

$(TYPEDEF)

Structures to contain the details of a solute or solvent to
store in the results of the MDDF calculation.

$(TYPEDFIELDS)

"""
struct SolSummary
    natoms::Int
    nmols::Int
    natomspermol::Int
end
SolSummary(s::Selection) = SolSummary(s.natoms, s.nmols, s.natomspermol)

#
# obs: voltar e deixar só a mutable struct
#
macro ResultFields_Start()
    ex = quote
        nbins::Int
        dbulk::Float64
        cutoff::Float64
        d::Vector{Float64} = zeros(nbins)

        # Data to compute the MDDF distribution and corresponding KB integral
        md_count::Vector{Float64} = zeros(nbins)
        md_count_random::Vector{Float64} = zeros(nbins)
        sum_md_count::Vector{Float64} = zeros(nbins)
        sum_md_count_random::Vector{Float64} = zeros(nbins)
        mddf::Vector{Float64} = zeros(nbins)
        kb::Vector{Float64} = zeros(nbins)

        # Properties of the solute and solvent selections
        autocorrelation::Bool
        solvent::SolSummary
        solute::SolSummary
    end
    esc(ex)
end

macro ResultFields_AtomsMatrix()
    ex = quote
        # Atomic contributions to the MDDFs
        solute_atom::Matrix{Float64} = zeros(nbins, solute.natomspermol)
        solvent_atom::Matrix{Float64} = zeros(nbins, solvent.natomspermol)
    end
    esc(ex)
end

macro ResultFields_AtomsVector()
    ex = quote
        # Atomic contributions to the MDDFs
        solute_atom::Array{Float64} = zeros(nbins, solute.natomspermol)
        solvent_atom::Array{Float64} = zeros(nbins, solvent.natomspermol)
    end
    esc(ex)
end

macro ResultFields_End()
    ex = quote
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
    esc(ex)
end

"""

$(TYPEDEF)

Structure to contain the results of the MDDF calculation.

$(TYPEDFIELDS)


"""
@with_kw_noshow struct Result
    @ResultFields_Start()
    @ResultFields_AtomsMatrix()
    @ResultFields_End()
end

# The mutable version is used for reading saved data, because some vectors
# need to be reshaped
@with_kw_noshow mutable struct MutableResult
    @ResultFields_Start()
    @ResultFields_AtomsVector()
    @ResultFields_End()
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
            error(
                " in MDDF options: The bulk volume is zero (dbulk must be smaller than cutoff). ",
            )
        end
        if (options.cutoff / options.binstep) % 1 > 1.e-5
            error("in MDDF options: cutoff must be a multiple of binstep.")
        end
        if ((options.cutoff - options.dbulk) / options.binstep) % 1 > 1.e-5
            error("in MDDF options: (cutoff-dbulk) must be a multiple of binstep. ")
        end
        cutoff = options.cutoff
    end
    nbins = setbin(cutoff, options.binstep) - 1

    if options.irefatom > trajectory.solvent.natoms
        error(
            "in MDDF options: Reference atom index",
            options.irefatom,
            " is greater than number of " *
            "                 atoms of the solvent molecule. ",
        )
    end

    # Set reference atom as the closest one to the center of coordinates of the molecule, as default
    if irefatom == -1
        if options.irefatom == -1
            nextframe!(trajectory)
            xfirst = trajectory.x_solvent[1:trajectory.solvent.natomspermol]
            cm = centerofcoordinates(xfirst)
            dmin, one, irefatom = minimumdistance(cm, xfirst)
            firstframe(trajectory)
        else
            irefatom = options.irefatom
        end
    end

    # Last frame to be considered
    if options.lastframe == -1
        lastframe_read = trajectory.nframes
    else
        lastframe_read = options.lastframe
    end

    # Actual number of frames that are read considering lastframe and stride
    nframes_read = round(Int, (lastframe_read - options.firstframe + 1) / options.stride)
    if nframes_read == 0
        error("Number of frames to read is zero. Check input parameters.")
    end

    # Return data structure built up

    return Result(
        options = options,
        nbins = nbins,
        dbulk = options.dbulk,
        cutoff = cutoff,
        irefatom = irefatom,
        lastframe_read = lastframe_read,
        nframes_read = nframes_read,
        autocorrelation = isequal(trajectory.solute.index,trajectory.solvent.index),
        solute = SolSummary(trajectory.solute),
        solvent = SolSummary(trajectory.solvent),
        files = [trajectory.filename],
        weights = [1.0],
    )

end

#
# What to show at the REPL
#
Base.show(io::IO, R::Result) = Base.show(io, overview(R))

"""

$(TYPEDEF)

Simple structure to contain the number of samples of each type of calculation to compute final results

$(TYPEDFIELDS)

"""
@with_kw struct Samples
    md::Float64
    random::Int
end

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

Compute the point in which the radius comprises half of the volume of the shell.

"""
function shellradius(i, step)
    rmin = (i - 1) * step
    return (0.5 * ((rmin + step)^3 + rmin^3))^(1 / 3)
end


"""
  sphereradiusfromshellvolume(volume,step)

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


"""
    finalresults!(R::Result, options::Options, trajectory::Trajectory, samples::Samples)

Function that computes the final results of all the data computed by averaging
according to the sampling of each type of data, and converts to common units.

Computes also the final distribution functions and KB integrals

This function modified the values contained in the R data structure

"""
function finalresults!(
    R::Result,
    options::Options,
    trajectory::Trajectory,
    samples::Samples,
)

    # Setup the distance vector
    for i = 1:R.nbins
        R.d[i] = shellradius(i, options.binstep)
    end

    # Counters
    @. R.md_count = R.md_count / (samples.md * R.nframes_read)
    @. R.solute_atom = R.solute_atom / (samples.md * R.nframes_read)
    @. R.solvent_atom = R.solvent_atom / (samples.md * R.nframes_read)
    @. R.md_count_random = R.md_count_random / (samples.random * R.nframes_read)
    @. R.rdf_count = R.rdf_count / (samples.md * R.nframes_read)
    @. R.rdf_count_random = R.rdf_count_random / (samples.random * R.nframes_read)

    # Volumes and Densities
    R.volume.total = R.volume.total / R.nframes_read
    R.density.solvent = R.density.solvent / R.nframes_read
    R.density.solute = R.density.solute / R.nframes_read

    R.volume.shell = R.volume.shell / R.nframes_read
    R.volume.domain = R.volume.domain / R.nframes_read
    R.volume.bulk = R.volume.bulk / R.nframes_read

    R.density.solvent_bulk = R.density.solvent_bulk / R.nframes_read

    #
    # Computing the distribution functions and KB integrals, from the MDDF
    # and from the RDF
    #

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
        end
        if ibin == 1
            R.sum_md_count[ibin] = R.md_count[ibin]
            R.sum_md_count_random[ibin] = R.md_count_random[ibin]
        else
            R.sum_md_count[ibin] = R.sum_md_count[ibin-1] + R.md_count[ibin]
            R.sum_md_count_random[ibin] =
                R.sum_md_count_random[ibin-1] + R.md_count_random[ibin]
        end
        R.kb[ibin] =
            units.Angs3tocm3permol *
            (1 / R.density.solvent_bulk) *
            (R.sum_md_count[ibin] - R.sum_md_count_random[ibin])

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

    return nothing
end

"""
    merge(r::Vector{ComplexMixutures.Result})

This function merges the results of MDDF calculations obtained by running the same
analysis on multiple trajectories, or multiple parts of the same trajectory. It returns
a Result structure of the same type, with all the functions and counters representing averages
of the set provided weighted by the number of frames read in each Result set.

"""
function Base.merge(r::Vector{Result})

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
        if (r[ir].cutoff - r[1].cutoff) > 1.e-8
            println(
                "ERROR: To merge Results, cutoff distance of the of the histograms of both sets must be the same.",
            )
        end
    end
    if error
        error(" Incompatible set of results to merge. ")
    end

    # List of files and weights
    nfiles = 0
    for ir = 1:nr
        nfiles += length(r[ir].files)
    end
    files = Vector{String}(undef, nfiles)
    weights = Vector{Float64}(undef, nfiles)

    # Final resuls
    R = Result(
        options = r[1].options,
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

    # Average results weighting the data considering the number of frames of each data set
    @. R.d = r[1].d
    ifile = 0
    for ir = 1:nr

        w = r[ir].nframes_read / nframes_read

        @. R.mddf += w * r[ir].mddf
        @. R.kb += w * r[ir].kb

        @. R.rdf += w * r[ir].rdf
        @. R.kb_rdf += w * r[ir].kb_rdf

        @. R.md_count += w * r[ir].md_count
        @. R.md_count_random += w * r[ir].md_count_random

        @. R.sum_md_count += w * r[ir].sum_md_count
        @. R.sum_md_count_random += w * r[ir].sum_md_count_random

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
            R.files[ifile] = r[ir].files[j]
            R.weights[ifile] = w * r[ir].weights[j]
        end

    end

    return R
end

#
# Functions to save the results to a file
#
StructTypes.StructType(::Type{SolSummary}) = StructTypes.Struct()
StructTypes.StructType(::Type{MutableResult}) = StructTypes.Struct()
StructTypes.StructType(::Type{Result}) = StructTypes.Struct()
StructTypes.StructType(::Type{Density}) = StructTypes.Struct()
StructTypes.StructType(::Type{Volume}) = StructTypes.Struct()
StructTypes.StructType(::Type{Options}) = StructTypes.Struct()

"""
    save(R::Result, filename::String)

Function to write the result data structure to a json file.

"""
function save(R::Result, filename::String)
    f = open(filename, "w")
    JSON3.write(f, R)
    close(f)
end


"""
  load(filename::String)

Function to load the json saved results file into the `Result` data structure.

"""
function load(filename::String)
    f = open(filename, "r")
    R = JSON3.read(f, MutableResult)
    # Need to reshape the solute and solvent atom contributions, because
    # the data is read in a single column
    solute_atom = reshape(R.solute_atom, R.nbins, :)
    solvent_atom = reshape(R.solvent_atom, R.nbins, :)
    R.solute_atom = solute_atom
    R.solvent_atom = solvent_atom
    return Result([getfield(R, field) for field in fieldnames(Result)]...)
end

# Format numbers depending on their size
format(x) = abs(x) < 999 ? @sprintf("%12.7f", x) : @sprintf("%12.5e", x)

import Base.write

"""
    write(R::ComplexMixtures.Result, filename::String, solute::Selection, solvent::Selection)

Function to write the final results to output files as simple tables that are human-readable and easy to analyze with other software

If the solute and solvent selections are provides, pass on the atom names.

"""
write(R::Result, filename::String, solute::Selection, solvent::Selection) =
    write(R, filename, solute_names = solute.names, solvent_names = solvent.names)


"""
    write(R::ComplexMixtures.Result, filename::String; 
          solute_names::Vector{String} = ["nothing"], 
          solvent_names::Vector{String} = ["nothing"])

Optional passing of atom names.

"""
function write(
    R::Result,
    filename::String;
    solute_names::Vector{String} = ["nothing"],
    solvent_names::Vector{String} = ["nothing"],
)

    # Names of output files containing atomic contibutions
    atom_contrib_solvent =
        FileOperations.remove_extension(filename) *
        "-ATOM_CONTRIB_SOLVENT." *
        FileOperations.file_extension(filename)
    atom_contrib_solute =
        FileOperations.remove_extension(filename) *
        "-ATOM_CONTRIB_SOLUTE." *
        FileOperations.file_extension(filename)


    #
    # GMD computed with minimum distance
    #

    # Conversion factor for volumes (as KB integrals), from A^3 to cm^3/mol
    mole = 6.022140857e23
    convert = mole / 1.e24

    output = open(filename, "w")
    println(output, @sprintf("#"))
    println(output, @sprintf("# Output of ComplexMixtures - MDDF"))
    println(output, @sprintf("# Trajectory files and weights:"))
    for i = 1:length(R.files)
        println(output, "#  $(R.files[i]) - w = $(R.weights[i])")
    end
    println(output, @sprintf("#"))
    println(
        output,
        @sprintf(
            "# Density of solvent in simulation box (sites/A^3): %15.8f",
            R.density.solvent
        )
    )
    println(
        output,
        @sprintf(
            "# Density of solvent in bulk (estimated) (sites/A^3): %15.8f",
            R.density.solvent_bulk
        )
    )
    println(
        output,
        @sprintf(
            "# Molar volume of solvent in simulation (cc/mol): %15.8f",
            convert / R.density.solvent
        )
    )
    println(
        output,
        @sprintf(
            "# Molar volume of solvent in bulk (estimated) (cc/mol): %15.8f",
            convert / R.density.solvent_bulk
        )
    )
    println(output, @sprintf("#"))
    println(output, @sprintf("# Number of atoms solute: %i9", size(R.solute_atom, 2)))
    println(
        output,
        @sprintf("# Number of atoms of the solvent: %i9", size(R.solvent_atom, 2))
    )
    println(output, @sprintf("#"))

    if R.options.usecutoff
        ibulk = setbin(R.options.dbulk, R.options.binstep)
    else
        ibulk = round(Int, R.nbins - 1 / R.options.binstep)
    end
    bulkerror = Statistics.mean(R.mddf[ibulk:R.nbins])
    sdbulkerror = Statistics.std(R.mddf[ibulk:R.nbins])
    println(output, "#")
    println(
        output,
        @sprintf(
            "# Average and standard deviation of bulk-gmd: %12.5f +/- %12.5f",
            bulkerror,
            sdbulkerror
        )
    )
    println(output, 
    """
    # COLUMNS CORRESPOND TO:
    #       1  Minimum distance to solute (dmin)
    #       2  GMD distribution (md count normalized by md count of random-solute distribution
    #       3  Kirwood-Buff integral (cc/mol) computed [(1/bulkdensity)*(col(6)-col(7))].
    #       4  Minimum distance site count for each dmin.
    #       5  Minimum distance site count for each dmin for random solute distribution.
    #       6  Cumulative number of molecules within dmin in the simulation.
    #       7  Cumulative number of molecules within dmin for random solute distribution.
    #       8  Volume of the shell of distance dmin and width binstep.
    #
    #   1-DISTANCE         2-GMD      3-KB INT    4-MD COUNT  5-COUNT RAND      6-SUM MD    7-SUM RAND   8-SHELL VOL
    """)

    for i = 1:R.nbins
        line = "  " * format(R.d[i])                                   #  1-DISTANCE
        line = line * "  " * format(R.mddf[i])                         #  2-GMD
        line = line * "  " * format(R.kb[i])                           #  3-KB INT
        line = line * "  " * format(R.md_count[i])                     #  4-MD COUNT
        line = line * "  " * format(R.md_count_random[i])              #  5-COUNT RAND
        line = line * "  " * format(R.sum_md_count[i])                 #  6-SUM MD
        line = line * "  " * format(R.sum_md_count_random[i])          #  7-SUM RAND
        line = line * "  " * format(R.volume.shell[i])                 #  8-SHELL VOL
        println(output, line)
    end
    close(output)

    # Writting gmd per atom contributions for the solvent

    output = open(atom_contrib_solvent, "w")
    println(output,
    """
    # Solvent atomic contributions to total MDDF.
    #
    # Trajectory files: $(R.files)
    #
    # Atoms: 
    """)
    for i = 1:size(R.solvent_atom, 2)
        if solvent_names[1] == "nothing"
            println(output, @sprintf("# %9i", i))
        else
            println(output, @sprintf("# %9i %5s", i, solvent_names[i]))
        end
    end
    println(output, "#")
    string = "#   DISTANCE     GMD TOTAL"
    for i = 1:size(R.solvent_atom, 2)
        string = string * @sprintf("  %12i", i)
    end
    println(output, string)
    for i = 1:R.nbins
        string = format(R.d[i])
        string = string * "  " * format(R.mddf[i])
        for j = 1:size(R.solvent_atom, 2)
            string = string * "  " * format(R.solvent_atom[i, j])
        end
        println(output, string)
    end
    close(output)

    # Writting gmd per atom contributions for the solute
    output = open(atom_contrib_solute, "w")
    println(output, 
    """
    #
    # Solute atomic contributions to total MDDF.
    #
    # Trajectory files: $(R.files)
    #
    # Atoms
    """)
    for i = 1:size(R.solute_atom, 2)
        if solute_names[1] == "nothing"
            println(output, @sprintf("# %9i", i))
        else
            println(output, @sprintf("# %9i %5s", i, solute_names[i]))
        end
    end
    println(output, "#")
    string = "#   DISTANCE      GMD TOTAL"
    for i = 1:size(R.solute_atom, 2)
        string = string * @sprintf("  %12i", i)
    end
    println(output, string)
    for i = 1:R.nbins
        string = format(R.d[i]) * "  " * format(R.mddf[i])
        for j = 1:size(R.solute_atom, 2)
            string = string * "  " * format(R.solute_atom[i, j])
        end
        println(output, string)
    end
    close(output)

    # Write final messages with names of output files and their content

    println()
    println(" OUTPUT FILES: ")
    println()
    println(" Wrote solvent atomic GMD contributions to file: ", atom_contrib_solvent)
    println(" Wrote solute atomic GMD contributions to file: ", atom_contrib_solute)
    println()
    println(" Wrote main output file: ", filename)
    println()

    return nothing
end

"""
  which_types(s::Selection, indexes::Vector{Int})

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
        if isel == nothing
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
  contrib(s::Selection, atom_contributions::Array{Float64}, selection)

Extract the contribution of a given atom type selection from the 
solute or solvent atomic contributions to the MDDF.

`s` here is the solute or solvent selection (type `ComplexMixtures.Selection`)
`atom_contributions` is the `R.solute_atom` or `R.solvent_atom` arrays of the `Result` structure,
and the last argument is the selection of atoms from the solute to be considered, given as a list 
of indexes, list of atom names, vector of `PDBTools.Atom`s, or a `PDBTools.Residue`. 

"""
function contrib(s::Selection, atom_contributions::Array{Float64}, indexes::Vector{Int})
    nbins = size(atom_contributions, 1)
    c = zeros(nbins)
    # If the selection is a single molecule, the indexes are anything
    if s.nmols == 1
        for it in indexes
            ind = findfirst(isequal(it),s.index)
            if isnothing(ind)
                error("Index $it of input list not found in selection indexes list.")
            end
            c += atom_contributions[:, ind]
        end
    # If more than one molecule, the index must correspond to an atom within one molecule
    else
        for it in indexes
            if it > s.natomspermol
                error(
                    "The index list contains atoms with indexes greater than the number of atoms of the molecule.",
                )
            end
            c += atom_contributions[:, it]
        end
    end
    return c
end

#
# If a list of atom names is provided
#
function contrib(s::Selection, atom_contributions::Array{Float64}, names::Vector{String})
    indexes = Vector{Int}(undef, 0)
    for name in names
        index = findall(isequal(name), s.names)
        if length(index) == 0
            error(" Atom in input list is not part of solvent (or solute): $name")
        end
        append!(indexes, index)
    end
    return contrib(s, atom_contributions, indexes)
end

#
# If a list of atoms of PDBTools.Atom is provided
#
function contrib(
    s::Selection,
    atom_contributions::Array{Float64},
    atoms::Vector{PDBTools.Atom};
    warning = true,
)
    (warning && s.nmols > 1) && warning_nmols_types()
    indexes = [atom.index for atom in atoms]
    # Check which types of atoms belong to this selection
    selected_types = which_types(s, indexes, warning = warning)
    return contrib(s, atom_contributions, selected_types)
end

#
# If a residue of type PDBTools.Residue is provided
#
function contrib(
    s::Selection,
    atom_contributions::Array{Float64},
    residue::Residue;
    warning = true,
)
    (warning && s.nmols > 1) && warning_nmols_types()
    indexes = collect(residue.range)
    # Check which types of atoms belong to this selection
    selected_types = which_types(s, indexes, warning = warning)
    return contrib(s, atom_contributions, selected_types)
end

function warning_nmols_types()
    println(
        "WARNING: there is more than one molecule in this selection. " *
        "Contributions are summed over all atoms of the same type.",
    )
end


"""
    sum!(R1::Result, R2::Result)

Sum the counts of two Results structures, adding the result to the first structure as in R1 = R1 + R2.

"""
function sum!(R1::Result, R2::Result)

    @. R1.md_count += R2.md_count
    @. R1.md_count_random += R2.md_count_random

    @. R1.solute_atom += R2.solute_atom
    @. R1.solvent_atom += R2.solvent_atom

    @. R1.rdf_count += R2.rdf_count
    @. R1.rdf_count_random += R2.rdf_count_random

    sum!(R1.density, R2.density)
    sum!(R1.volume, R2.volume)

    return nothing
end

function sum!(D1::Density, D2::Density)
    D1.solute += D2.solute
    D1.solvent += D2.solvent
    D1.solvent_bulk += D2.solvent_bulk
    return nothing
end

function sum!(V1::Volume, V2::Volume)
    V1.total += V2.total
    V1.bulk += V2.bulk
    V1.domain += V2.domain
    @. V1.shell += V2.shell
    return nothing
end

"""
    title(R::Result, solute::Selection, solvent::Selection)
    title(R::Result, solute::Selection, solvent::Selection, nspawn::Int)

$(INTERNAL)

Print some information about the run.

"""
function title(R::Result, solute::Selection, solvent::Selection)
    print("""
          $(bars)
          Starting MDDF calculation:
          $(R.nframes_read) frames will be considered.
          Solute: $(atoms_str(solute.natoms)) belonging to $(mol_str(solute.nmols)).
          Solvent: $(atoms_str(solvent.natoms)) belonging to $(mol_str(solvent.nmols))
          """)
end
function title(R::Result, solute::Selection, solvent::Selection, nspawn::Int)
    print(""" 
          $(bars)
          Starting MDDF calculation in parallel:
          $(R.nframes_read) frames will be considered.
          Number of calculation threads: $(nspawn)
          Solute: $(atoms_str(solute.natoms)) belonging to $(mol_str(solute.nmols)).
          Solvent: $(atoms_str(solvent.natoms)) belonging to $(mol_str(solvent.nmols)).
          """)
end
