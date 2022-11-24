"""

$(TYPEDEF)

Structure to contain PDB trajectories. Frames must be separated by "END", and with periodic cell sizes
in the "CRYST1" field.

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

    # This vector must be filled up with the size vectors of the periodic cell
    sides::Vector{T}

    # Solute and solvent data
    solute::Selection
    solvent::Selection

    # Coordinates of the solute and solvent atoms in a frame (natoms,3) for each array:
    x_solute::Vector{T}  # solute.natoms vectors of length 3
    x_solvent::Vector{T} # solvent.natoms vectors of length 3

end

"""
    PDBTraj(pdbfile::String, solute::Selection, solvent::Selection;T::Type = SVector{3,Float64})

Function open will set up the IO stream of the trajectory, fill up the number of frames field and additional parameters if required 

"""
function PDBTraj(
    pdbfile::String,
    solute::Selection,
    solvent::Selection;
    T::Type = SVector{3,Float64},
)

    st = open(pdbfile, "r")

    # Get the number of frames and number of atoms, check for adequacy of the format
    nframes = 0
    natoms = 0
    for line in eachline(st)
        isempty(line) && continue
        line_data = split(line)
        if line_data[1] == "END"
            nframes = nframes + 1
        end
        if nframes == 0 && line_data[1] in ("ATOM","HETATM")
            natoms = natoms + 1
        end
    end
    close(st)

    # Some error messages
    if nframes == 0 
        error("Could not read any frame from PDB file. Each frame must end with the END specifier")
    end
    if natoms == 0
        error("Could not read any ATOM from the trajectory file.")
    end

    # Fill-up the sides vector of the trajectory. We assume here
    # that the sides are stored for each frame in the "CRYST1" field, for each frame.
    # Here, we exemplify the option to read all sides at once and store them in a 
    # vector. Alternatively, the sides could be read for each frame
    # independently within the "nextframe!" function, and saved as a sides(3) vector.  
    # The function "getsides", below, must be adapted accordingly to return the correct
    # sides of the periodic box in each frame.
    sides = zeros(T, nframes)
    st = open(pdbfile, "r")
    iframe = 0
    for line in eachline(st)
        isempty(line) && continue
        s = split(line)
        if s[1] == "CRYST1"
            iframe = iframe + 1
            sides[iframe] = T(parse.(Float64, @view(s[2:4])))
        end
    end

    # Setup the struct that contains the stream in the trajectory type
    stream = Stream(st)

    # Return the stream closed, it is opened and closed within the mddf routine
    close(st)

    return PDBTraj(
        pdbfile, # trajectory file name 
        stream,
        nframes,
        sides, # array containing box sides
        solute,
        solvent,
        zeros(T, solute.natoms),
        zeros(T, solvent.natoms),
    )
end

function Base.show(io::IO, trajectory::PDBTraj)
    print(io,""" 
          Trajectory in PDB format with:
              $(trajectory.nframes) frames.
              Solute contains $(trajectory.solute.natoms) atoms.
              Solvent contains $(trajectory.solvent.natoms) atoms.
          """)
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
    st = stream(trajectory)
    iatom = 0
    record = readline(st)
    i_solute = 0
    i_solvent = 0
    while ((0 < length(record) < 3) || record[1:3] != "END")
        if length(record) >= 6
            if record[1:4] == "ATOM" || record[1:6] == "HETATM"
                iatom = iatom + 1
                x = parse(Float64, record[31:38])
                y = parse(Float64, record[39:46])
                z = parse(Float64, record[47:54])
                if i_solute < trajectory.solute.natoms && iatom == trajectory.solute.index[i_solute+1]
                    i_solute += 1
                    trajectory.x_solute[i_solute] = T(x,y,z)
                end
                if i_solvent < trajectory.solvent.natoms && iatom == trajectory.solvent.index[i_solvent+1]
                    i_solvent += 1
                    trajectory.x_solvent[i_solvent] = T(x,y,z)
                end
            end
        end
        record = readline(st)
    end
end

# Function that returns the sides of the periodic box given the data structure. In this
# case, returns the 3-element vector corresponding to the box sides of the given frame
function getsides(trajectory::PDBTraj, iframe)
    # Sides is expected to be an array that contains the sides for each frame, and we return the
    # vector containing the sides of the current fraem
    return trajectory.sides[iframe]
end

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
