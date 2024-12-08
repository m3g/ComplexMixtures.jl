#
# Structure to contain trajectories read by the Chemfiles package. Many formats should
# be available through this interface, including the NamdDCD which are provided independently
# as an example. 
#
import Chemfiles

# Add this format to list of available formats 
push!(trajectory_formats, "Chemfiles")

#"""
#
#$(TYPEDEF)
#
#Structure to contain a trajectory as read by Chemfiles.jl
#
#$(TYPEDFIELDS)
#
#"""
struct ChemFile{T<:AbstractVector} <: Trajectory

    #
    # Mandatory data for things to work
    #
    filename::String
    format::AbstractString
    stream::Stream{<:Chemfiles.Trajectory} # mutable such that we can close it and open it again
    nframes::Int64

    # Solute and solvent data
    solute::AtomSelection
    solvent::AtomSelection

    # Coordinates of the solute and solvent atoms in a frame (natoms,3) for each array:
    x_solute::Vector{T}  # natoms(solute) vectors of length 3 (preferentially static vectors)
    x_solvent::Vector{T} # natoms(solvent) vectors of lenght 3 (preferentially static vectors)

    # unitcell
    unitcell::MMatrix{3,3,Float64,9}

    #
    # Additional data required for input/output functions
    #
    natoms::Int64

end

#"""
#    ChemFile(filename::String, solute::AtomSelection, solvent::AtomSelection;format="" , T::Type = SVector{3,Float64})
#
#Function open will set up the IO stream of the trajectory, fill up the number of frames field and additional parameters if required.
#
#"""
function ChemFile(
    filename::String,
    solute::AtomSelection,
    solvent::AtomSelection;
    format="",
    T::Type=SVector{3,Float64},
)

    st = redirect_stdout(() -> Chemfiles.Trajectory(filename, 'r', format), devnull)

    # Get the number of frames (the output of Chemfiles comes in UInt64 format, which is converted
    # to Int using (UInt % Int)
    nframes = Chemfiles.length(st) % Int

    # Read the first frame to get the number of atoms
    frame = Chemfiles.read(st)
    natoms = Chemfiles.size(frame) % Int
    unitcell = transpose(MMatrix{3,3}(Chemfiles.matrix(Chemfiles.UnitCell(frame))))

    # Initialize the stream struct of the Trajectory
    stream = Stream(st)

    # Return the stream closed, it is opened and closed within the mddf routine
    Chemfiles.close(st)

    return ChemFile(
        filename, # trajectory file name 
        format, # trajectory format, is provided by the user
        stream,
        nframes,
        solute,
        solvent,
        zeros(T, ComplexMixtures.natoms(solute)),
        zeros(T, ComplexMixtures.natoms(solvent)),
        unitcell,
        natoms, # Total number of atoms
    )
end

function Base.show(io::IO, trajectory::ChemFile)
    print(io, strip("""
           Trajectory read by Chemfiles with:
               $(trajectory.nframes) frames.
               $(trajectory.natoms) atoms.
               Unit cell in current frame: $(convert_unitcell(getunitcell(trajectory)))
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
function nextframe!(trajectory::ChemFile{T}) where {T}

    st = stream(trajectory)

    frame = Chemfiles.read(st)
    positions = Chemfiles.positions(frame)
    trajectory.unitcell .= transpose(SMatrix{3,3}(Chemfiles.matrix(Chemfiles.UnitCell(frame))))

    # Save coordinates of solute and solvent in trajectory arrays (of course this could be avoided,
    # but the code in general is more clear aftwerwards by doing this)
    for i in eachindex(trajectory.x_solute)
        trajectory.x_solute[i] = T(
            positions[1, trajectory.solute.indices[i]],
            positions[2, trajectory.solute.indices[i]],
            positions[3, trajectory.solute.indices[i]],
        )
    end
    for i in eachindex(trajectory.x_solvent)
        trajectory.x_solvent[i] = T(
            positions[1, trajectory.solvent.indices[i]],
            positions[2, trajectory.solvent.indices[i]],
            positions[3, trajectory.solvent.indices[i]],
        )
    end

    return trajectory
end

# Returns the unitcell of the current frame, as a 3x3 matrix (static in this case)
getunitcell(trajectory::ChemFile) = trajectory.unitcell

#
# Function that closes the IO Stream of the trajectory
#
closetraj!(trajectory::ChemFile) = Chemfiles.close(stream(trajectory))

#
# Function to open the trajectory stream
#
function opentraj!(trajectory::ChemFile)
    st = redirect_stdout(
        () -> Chemfiles.Trajectory(trajectory.filename, 'r', trajectory.format),
        devnull,
    )
    set_stream!(trajectory, st)
end

#
# Function that returns the trajectory in position to read the first frame
#
function firstframe!(trajectory::ChemFile)
    closetraj!(trajectory)
    opentraj!(trajectory)
end

@testitem "getunitcell" begin
    using ComplexMixtures
    using ComplexMixtures.Testing
    using PDBTools
    using StaticArrays

    atoms = readPDB(Testing.pdbfile)
    options = Options(stride=4, seed=321, StableRNG=true, nthreads=1, silent=true)
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    traj = Trajectory("$(Testing.data_dir)/NAMD/trajectory.dcd", protein, tmao)
    ComplexMixtures.opentraj!(traj)
    ComplexMixtures.firstframe!(traj)
    ComplexMixtures.nextframe!(traj)
    @test ComplexMixtures.convert_unitcell(ComplexMixtures.getunitcell(traj)) ≈
          SVector(84.42188262939453, 84.42188262939453, 84.42188262939453)
    ComplexMixtures.closetraj!(traj)
end
