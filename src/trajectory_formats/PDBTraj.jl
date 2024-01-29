"""

$(TYPEDEF)

Structure to contain PDB trajectories. Frames must be separated by "END", and with periodic cell sizes
in the "CRYST1" field, for each frame.

This structure and functions can be used as a template to implement the reading of other trajectory formats. 

$(TYPEDFIELDS)

"""
struct PDBTraj{T<:AbstractVector} <: Trajectory

    #
    # Mandatory data for things to work
    #
    filename::String
    stream::Stream{<:IOStream} # The type of IOStream must be the result of open(filename, "r") for the given type
    nframes::Int64

    # unitcell in the current frame
    unitcell::MMatrix{3,3,Float64,9}

    # Solute and solvent data
    solute::AtomSelection
    solvent::AtomSelection

    # Coordinates of the solute and solvent atoms in a frame (natoms,3) for each array:
    x_solute::Vector{T}  # solute number of atoms vectors of length 3
    x_solvent::Vector{T} # solvent number of atoms vectors of length 3

end

"""
    PDBTraj(pdbfile::String, solute::AtomSelection, solvent::AtomSelection;T::Type = SVector{3,Float64})

Function open will set up the IO stream of the trajectory, fill up the number of frames field and additional parameters if required 

"""
function PDBTraj(
    pdbfile::String,
    solute::AtomSelection,
    solvent::AtomSelection;
    T::Type = SVector{3,Float64},
)

    st = open(pdbfile, "r")

    # Get the number of frames and number of atoms, check for adequacy of the format
    nframes = 0
    natoms = 0
    unitcell = zeros(MMatrix{3,3,Float64,9})
    for line in eachline(st)
        isempty(line) && continue
        line_data = split(line)
        if line_data[1] == "CRYST1"
            unitcell .= pdb_readunitcell(line_data) 
        end
        if line_data[1] == "END"
            nframes = nframes + 1
        end
        if nframes == 0 && line_data[1] in ("ATOM", "HETATM")
            natoms = natoms + 1
        end
    end
    close(st)

    # Some error messages
    if all(unitcell .== 0)
        error("Could not read unit cell from PDB file. Each frame must contain a CRYST1 field.")
    end
    if nframes == 0
        error("Could not read any frame from PDB file. Each frame must end with the END specifier")
    end
    if natoms == 0
        error("Could not read any ATOM from the trajectory file.")
    end

    # Setup the struct that contains the stream in the trajectory type
    stream = Stream(st)

    # Return the stream closed, it is opened and closed within the mddf routine
    close(st)

    return PDBTraj(
        pdbfile, # trajectory file name 
        stream,
        nframes,
        unitcell, # unitcell in the current frame
        solute,
        solvent,
        zeros(T, ComplexMixtures.natoms(solute)),
        zeros(T, ComplexMixtures.natoms(solvent)),
    )
end

function Base.show(io::IO, trajectory::PDBTraj)
    (; solute, solvent) = trajectory
    print(io,strip(""" 
          Trajectory in PDB format with:
              $(trajectory.nframes) frames.
              Solute contains $(ComplexMixtures.natoms(solute)) atoms.
              Solvent contains $(ComplexMixtures.natoms(solvent)) atoms.
              Unit cell in current frame: $(print_unitcell(trajectory))
          """))
end

#
# Function that reads the coordinates of the solute and solvent atoms from
# the next frame of the trajectory file 
#
# The function modifies sides, x_solute and x_solvent within the trajectory structure.
# Having these vectors inside the trajectory structure avoids having to allocate
# them everytime a new frame is read
#
function nextframe!(trajectory::PDBTraj{T}) where {T<:AbstractVector}
    (; solute, solvent) = trajectory
    st = stream(trajectory)
    iatom = 0
    record = readline(st)
    i_solute = 0
    i_solvent = 0
    while ((0 < length(record) < 3) || record[1:3] != "END")
        if length(record) >= 6
            # Read unit cell
            if record[1:6] == "CRYST1"
                trajectory.unitcell .= pdb_readunitcell(split(record)) 
            end
            # Read atom coordinates into the solute and solvent arrays
            if record[1:4] == "ATOM" || record[1:6] == "HETATM"
                iatom = iatom + 1
                x = parse(Float64, record[31:38])
                y = parse(Float64, record[39:46])
                z = parse(Float64, record[47:54])
                if i_solute < ComplexMixtures.natoms(solute) &&
                   iatom == trajectory.solute.indices[i_solute+1]
                    i_solute += 1
                    trajectory.x_solute[i_solute] = T(x, y, z)
                end
                if i_solvent < ComplexMixtures.natoms(solvent) &&
                   iatom == trajectory.solvent.indices[i_solvent+1]
                    i_solvent += 1
                    trajectory.x_solvent[i_solvent] = T(x, y, z)
                end
            end
        end
        record = readline(st)
    end
end

#
# Function that returns the unitcell of the current trajctory frame
#
getunitcell(trajectory::PDBTraj) = trajectory.unitcell

#
# Function that opens the trajectory stream
#
opentraj!(trajectory::PDBTraj) = set_stream!(trajectory, open(trajectory.filename, "r"))

#
# Function that closes the IO Stream of the trajectory
#
closetraj!(trajectory::PDBTraj) = close(stream(trajectory))

#
# Function that returns the trajectory in position to read the first frame
#
firstframe!(trajectory::PDBTraj) = seekstart(stream(trajectory))

#
# Function to read a unit cell from the CRYST1 field of the PDB file
# 
# COLUMNS       DATA  TYPE    FIELD          DEFINITION
# -------------------------------------------------------------
#  1 -  6       Record name   "CRYST1"
#  7 - 15       Real(9.3)     a              a (Angstroms).
# 16 - 24       Real(9.3)     b              b (Angstroms).
# 25 - 33       Real(9.3)     c              c (Angstroms).
# 34 - 40       Real(7.2)     alpha          alpha (degrees).
# 41 - 47       Real(7.2)     beta           beta (degrees).
# 48 - 54       Real(7.2)     gamma          gamma (degrees).
# 56 - 66       LString       sGroup         Space  group.
# 67 - 70       Integer       z              Z value.
function pdb_readunitcell(data::AbstractVector{<:AbstractString})
    a = parse(Float64, data[2])
    b = parse(Float64, data[3])
    c = parse(Float64, data[4])
    alpha = parse(Float64, data[5])
    beta = parse(Float64, data[6])
    gamma = parse(Float64, data[7])
    return transpose(SMatrix{3,3}(Chemfiles.matrix(Chemfiles.UnitCell([a, b, c], [alpha, beta, gamma]))))
end
