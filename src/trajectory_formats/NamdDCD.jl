#
# Structure to contain DCD trajectories produces with Namd. 
#
"""

$(TYPEDEF)

Structure to contain the data of a trajectory in NAMD/DCD format.

$(TYPEDFIELDS)

"""
struct NamdDCD{T<:AbstractVector} <: Trajectory

    #
    # Mandatory data for things to work
    #
    filename::String
    stream::Stream{<:FortranFile} # special type of stream required for reading DCD files
    nframes::Int64

    # Data structures of the solute and solvent 
    solute::AtomSelection
    solvent::AtomSelection

    # Coordinates of the solute and solvent atoms in a frame (3,natoms) for each array:
    x_solute::Vector{T}
    x_solvent::Vector{T}

    #
    # Additional properties that might be required for implementing IO (not necessary for every
    # input/output format)
    #
    lastatom::Int64 # The last atom to be read from each line

    # Auxiliary vectors to read the coordinates without having the allocate/deallocate every time
    unitcell_read::Vector{Float64}
    x_read::Vector{Float32}
    y_read::Vector{Float32}
    z_read::Vector{Float32}

end

"""
    NamdDCD(filename::String, solute::AtomSelection, solvent::AtomSelection;T::Type = SVector{3,Float64})

This function initializes the structure above, returning the data and the vectors with appropriate lengths.

"""
function NamdDCD(
    filename::String,
    solute::AtomSelection,
    solvent::AtomSelection;
    T::Type = SVector{3,Float64},
)

    st = FortranFile(filename)

    # Read header
    IntVec = Vector{Int32}(undef, 17)
    _, _, IntVec[1:8], _, IntVec[9:17] = read(st, FString{4}, Int32, (Int32, 8), Float64, (Int32, 9))
    _, _ = read(st, Int32, FString{80})
    read_natoms = read(st, Int32)

    # Check if dcd file contains unit cell information
    unitcell_in_dcd = false
    try
        read(st, [Float64 for i = 1:7])
    catch err
        unitcell_in_dcd = true
    end
    if !unitcell_in_dcd
        error("DCD file does not contain unit cell information.")
    end

    # Get number of frames
    firstframe!(st)

    # Read the unitcell in the first frame
    unitcell_read = zeros(Float64, 6)
    read(st, unitcell_read)
    firstframe!(st)

    nframes = getnframes(st)
    lastatom = max(maximum(solute.index), maximum(solvent.index))

    # setup the trajectory struct that contains the stream
    stream = Stream(st)

    # Return the stream closed, it is opened and closed within the mddf routine
    FortranFiles.close(st)

    return NamdDCD(
        filename,
        stream,
        nframes,
        solute,
        solvent,
        zeros(T, solute.natoms), # solute atom coordinates
        zeros(T, solvent.natoms), # solvent atom coordinates
        lastatom,
        unitcell_read, # auxiliary vector to read unitcell data
        zeros(Float32, lastatom), # auxiliary x
        zeros(Float32, lastatom), # auxiliary y
        zeros(Float32, lastatom),  # auxiliary z
    )

end

function Base.show(io::IO, trajectory::NamdDCD)
    print(io, strip(""" 
          Trajectory in NamdDCD format containing:
              $(trajectory.nframes) frames.
              Unit cell in current frame: $(print_unitcell(trajectory))
          """))
end

#
# Function that opens the trajectory stream
#
opentraj!(trajectory::NamdDCD) = set_stream!(trajectory, FortranFile(trajectory.filename))

#
# Function that closes the IO Stream of the trajectory
#
closetraj!(trajectory::NamdDCD) = FortranFiles.close(stream(trajectory))

#
# Function that reads the coordinates of the solute and solvent atoms from
# the next frame of the trajectory file 
#
# The function modifies unitcell, x_solute and x_solvent within the trajectory structure.
# Having these vectors inside the trajectory structure avoids having to allocate
# them everytime a new frame is read
#
function nextframe!(trajectory::NamdDCD{T}) where {T}

    st = stream(trajectory)

    read(st, trajectory.unitcell_read)

    # Read the coordinates  
    read(st, trajectory.x_read)
    read(st, trajectory.y_read)
    read(st, trajectory.z_read)

    # Save coordinates of solute and solvent in trajectory arrays
    for i = 1:trajectory.solute.natoms
        trajectory.x_solute[i] = T(
            trajectory.x_read[trajectory.solute.index[i]],
            trajectory.y_read[trajectory.solute.index[i]],
            trajectory.z_read[trajectory.solute.index[i]],
        )
    end
    for i = 1:trajectory.solvent.natoms
        trajectory.x_solvent[i] = T(
            trajectory.x_read[trajectory.solvent.index[i]],
            trajectory.y_read[trajectory.solvent.index[i]],
            trajectory.z_read[trajectory.solvent.index[i]],
        )
    end

    return nothing
end

#
# Function that returns the unitcell given the way that the PBC information is 
# stored in the DCD trajectory structure.
#
function getunitcell(trajectory::NamdDCD)
    A = trajectory.unitcell_read[1]
    γ = trajectory.unitcell_read[2]
    B = trajectory.unitcell_read[3]
    β = trajectory.unitcell_read[4]
    α = trajectory.unitcell_read[5]
    C = trajectory.unitcell_read[6]
    if any(==(0.0), (α, β, γ))
        α = 90.0
        β = 90.0
        γ = 90.0
    end
    return transpose(SMatrix{3,3}(Chemfiles.matrix(Chemfiles.UnitCell([A,B,C], [α, β, γ]))))  
end

#
# Leave DCD file in position to read the first frame: DCD files have a header
#
function firstframe!(st::FortranFile)
    # rewind
    rewind(st)
    # skip header
    read(st)
    read(st)
    read(st)
end
firstframe!(trajectory::NamdDCD) = firstframe!(stream(trajectory))

#
# Auxiliary functions
#

#
# Sometimes the DCD files contains a wrong number of frames in the header, so to
# get the actual number of frames, it is better to read it
#
function getnframes(st::FortranFile)
    firstframe!(st)
    nframes = 0
    while true
        try
            read(st, Float64) # pbc data
            read(st, Float32) # x
            read(st, Float32) # y
            read(st, Float32) # z
            nframes = nframes + 1
        catch
            firstframe!(st)
            return nframes
        end
    end
end
getnframes(traj::NamdDCD) = getnframes(stream(traj))
