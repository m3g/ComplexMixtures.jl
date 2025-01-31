"""
    Trajectory(filename::String, solute::AtomSelection, solvent::AtomSelection; format::String = "", chemfiles = false)

Trajectory constructor data type. 

Defaults to reading with the Chemfiles infrastructure, except for DCD and PDB trajectory files, if the "PDBTraj" option is provided.

See memory issue (https://github.com/chemfiles/Chemfiles.jl/issues/44)

"""
abstract type Trajectory end

# Mutable struct to contain the stream of the trajectory: mutable such that
# we can open and close the trajectories
mutable struct Stream{T}
    st::T
end

# Get property with simplified syntax
stream(traj::Trajectory) = traj.stream.st
set_stream!(traj::Trajectory, st) = traj.stream.st = st

# Trajectory formats
const trajectory_formats = String[]

# Include specific predefined trajectory types
include("./trajectory_formats/ChemFiles.jl")
include("./trajectory_formats/NamdDCD.jl")
include("./trajectory_formats/PDBTraj.jl")

function Trajectory(
    filename::String,
    solute::AtomSelection,
    solvent::AtomSelection;
    format::String="",
    chemfiles=false,
)
    if !isempty(format) && !(format in trajectory_formats)
        throw(ArgumentError("""\n
            Trajectory format not properly set: $format 
            Available trajectory formats: $(join(trajectory_formats, ", "))

        """))
    end
    filename = expanduser(filename) # expand tilde on Unix systems, to username
    file_extension = filename[findlast(==('.'), filename)+1:end]
    trajectory = if !chemfiles && (format == "dcd" || file_extension == "dcd")
        NamdDCD(filename, solute, solvent)
    elseif !chemfiles && (format == "PDBTraj" || file_extension =="pdb")
        PDBTraj(filename, solute, solvent)
    else
        ChemFile(filename, solute, solvent, format=format)
    end
    return trajectory
end

# If only one selection is provided, assume that the solute and the solvent are the same
Trajectory(filename::String, solvent::AtomSelection; format::String="", chemfiles=false) =
    Trajectory(filename, solvent, solvent; format, chemfiles)

#=
    convert_unitcell(unitcell::Union{SVector{3}, SMatrix{3,3}})

Function to return the unit cell as a vector or matrix, depending on if 
the cell is diagonal or not, up to a relative precision of 1e-10 by default.

=#
function convert_unitcell(unitcell::AbstractMatrix; tol=1e-10)
    size(unitcell) == (3, 3) || error("Unit cell must be a 3x3 matrix.")
    s = minimum(diag(unitcell))
    is_diag = all(unitcell[i, j] < tol * s for i in 1:3, j in 1:3 if i != j)
    return is_diag ? SVector{3}(diag(unitcell)) : SMatrix{3,3}(unitcell)
end

# Version to ensure type stability when we know the type of the unit cell
convert_unitcell(::SVector, unitcell::AbstractMatrix) = SVector{3}(unitcell[1, 1], unitcell[2, 2], unitcell[3, 3])
convert_unitcell(::SMatrix, unitcell::AbstractMatrix) = SMatrix{3,3}(unitcell)

@testitem "convert_unitcell" setup=[AllocTest] begin
    using ComplexMixtures: convert_unitcell
    using StaticArrays: SVector, SMatrix
    using BenchmarkTools: @benchmark
    using .AllocTest: Allocs

    m = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test convert_unitcell(m) isa SVector
    s = SVector(1.0, 1.0, 1.0)
    @test convert_unitcell(s, m) isa SVector
    @test convert_unitcell(s, m) isa SVector
    b = @benchmark convert_unitcell($s, $m) samples = 1 evals = 1
    @test b.allocs ==  Allocs(0)
    m = [1.0 1.0 0.0; 0.0 1.0 0.0; 0.0 0.0 2.0]
    @test convert_unitcell(m) isa SMatrix
    m2 = SMatrix{3,3}(m)
    @test convert_unitcell(m2, m) isa SMatrix
    b = @benchmark convert_unitcell($m2, $m) samples = 1 evals = 1
    @test b.allocs == Allocs(0)
end

function print_unitcell(trajectory)
    unitcell = convert_unitcell(getunitcell(trajectory))
    if unitcell isa SVector
        return "[ $(@sprintf("%6.2f", unitcell[1])) 0 0;" *
               " 0 $(@sprintf("%6.2f", unitcell[2])) 0;" *
               " 0 0 $(@sprintf("%6.2f", unitcell[3])) ]"
    else
        return "[ $(@sprintf("%6.2f %6.2f %6.2f", unitcell[1,1], unitcell[1,2], unitcell[1,3]));" *
               " $(@sprintf("%6.2f %6.2f %6.2f", unitcell[2,1], unitcell[2,2], unitcell[2,3]));" *
               " $(@sprintf("%6.2f %6.2f %6.2f", unitcell[3,1], unitcell[3,2], unitcell[3,3])) ]"
    end
end

@testitem "Trajectory" begin
    using ComplexMixtures
    using ComplexMixtures.Testing
    using PDBTools
    using StaticArrays

    atoms = readPDB(Testing.pdbfile)
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)

    # Throw if trajectory type is not available
    @test_throws ArgumentError Trajectory("$(Testing.data_dir)/PDB/trajectory.pdb", protein, tmao; format="XXX")

    # NAMD DCD file
    traj = Trajectory("$(Testing.data_dir)/NAMD/trajectory.dcd", protein, tmao)
    @test traj.nframes == 20
    @test traj.lastatom == 4012
    @test ComplexMixtures.convert_unitcell(ComplexMixtures.getunitcell(traj)) ≈ SVector(84.42188262939453, 84.42188262939453, 84.42188262939453)

    # PDB file
    traj = Trajectory(
        "$(Testing.data_dir)/PDB/trajectory.pdb",
        protein,
        tmao,
        format="PDBTraj",
    )
    @test ComplexMixtures.natoms(traj.solute) == 1463
    @test ComplexMixtures.natoms(traj.solvent) == 2534
    @test ComplexMixtures.convert_unitcell(ComplexMixtures.getunitcell(traj)) ≈ SVector(84.47962951660156, 84.47962951660156, 84.47962951660156)
    # Test determining format from file extension
    traj = Trajectory("$(Testing.data_dir)/PDB/trajectory.pdb", protein, tmao)
    @test ComplexMixtures.natoms(traj.solute) == 1463

    # Chemfiles with NAMD
    traj = Trajectory(
        "$(Testing.data_dir)/NAMD/trajectory.dcd",
        protein,
        tmao,
        chemfiles=true,
    )
    @test traj.nframes == 20
    @test ComplexMixtures.convert_unitcell(ComplexMixtures.getunitcell(traj)) ≈ SVector(84.42188262939453, 84.42188262939453, 84.42188262939453)
    @test ComplexMixtures.natoms(traj.solute) == 1463
    @test ComplexMixtures.natoms(traj.solvent) == 2534

    # Chemfiles with Gromacs
    atoms = readPDB("$(Testing.data_dir)/Gromacs/system.pdb")
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    emi = AtomSelection(select(atoms, "resname EMI"), natomspermol=20)
    traj = Trajectory("$(Testing.data_dir)/Gromacs/trajectory.xtc", protein, emi)
    @test traj.nframes == 26
    @test ComplexMixtures.convert_unitcell(ComplexMixtures.getunitcell(traj)) ≈ SVector(95.11481285095215, 95.11481285095215, 95.13440132141113)
    @test ComplexMixtures.natoms(traj.solute) == 1231
    @test ComplexMixtures.natoms(traj.solvent) == 5080

end

#
# This structure and function are used to retrieve metadata from the trajectory
# that is used to initialize the Result and ParticleSystem structures, without
# having to open the trajectory again. This allows the parallel initialization of
# these data structures. 
#
@kwdef struct TrajectoryMetaData{UC}
    irefatom::Int
    lastframe_read::Int
    nframes_read::Int
    n_groups_solute::Int
    n_groups_solvent::Int
    unitcell::UC
end

function TrajectoryMetaData(trajectory::Trajectory, options::Options)

    if options.irefatom > trajectory.solvent.natomspermol
        throw(ArgumentError("in MDDF options: Reference atom index $(options.irefatom) is greater than number of atoms of the solvent molecule. "))
    end
    if options.lastframe > trajectory.nframes
        throw(ArgumentError("in MDDF options: lastframe is greater than trajectory.nframes. "))
    end

    # Open trajectory to read some data
    opentraj!(trajectory)
    firstframe!(trajectory)

    # Get unitcell from the trajectory: returns vector or matrix depending on the data
    unitcell = convert_unitcell(getunitcell(trajectory))

    # Set reference atom as the closest one to the center of coordinates of the molecule, as default
    if options.irefatom == -1
        nextframe!(trajectory)
        first_mol = viewmol(1, trajectory.x_solvent, trajectory.solvent)
        cm = mean(first_mol)
        irefatom = last(findmin(at -> norm(at - cm), first_mol))
    else
        irefatom = options.irefatom
    end

    # Last frame to be considered
    if options.lastframe == -1
        lastframe_read = trajectory.nframes
    else
        lastframe_read = options.lastframe
    end

    # Actual number of frames that are read considering lastframe and stride
    nframes_read = length(options.firstframe:options.stride:lastframe_read)

    # Close trajecotory
    closetraj!(trajectory)

    # Initialize the arrays that contain groups counts, depending on wheter
    # groups were defined or not in the input Options
    n_groups_solute = if !trajectory.solute.custom_groups
        trajectory.solute.natomspermol
    else
        length(trajectory.solute.group_atom_indices)
    end
    n_groups_solvent = if !trajectory.solvent.custom_groups
        trajectory.solvent.natomspermol
    else
        length(trajectory.solvent.group_atom_indices)
    end

    return TrajectoryMetaData(
        irefatom=irefatom,
        lastframe_read=lastframe_read,
        nframes_read=nframes_read,
        n_groups_solute=n_groups_solute,
        n_groups_solvent=n_groups_solvent,
        unitcell=unitcell,
    )
end

@testitem "TrajectoryMetaData" begin
    using ComplexMixtures
    using ComplexMixtures.Testing
    using PDBTools
    atoms = readPDB(Testing.pdbfile)
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    traj = Trajectory("$(Testing.data_dir)/NAMD/trajectory.dcd", protein, tmao)

    options = Options()
    tmeta = ComplexMixtures.TrajectoryMetaData(traj, options)
    @test tmeta.irefatom == 1
    @test tmeta.lastframe_read == 20
    @test tmeta.nframes_read == 20
    @test tmeta.n_groups_solute == 1463
    @test tmeta.n_groups_solvent == 14
    @test tmeta.unitcell ≈ [84.42188262939453, 84.42188262939453, 84.42188262939453]

    options = Options(; irefatom=2, lastframe=10, stride=2)
    tmeta = ComplexMixtures.TrajectoryMetaData(traj, options)
    @test tmeta.irefatom == 2
    @test tmeta.lastframe_read == 10
    @test tmeta.nframes_read == 5
    @test tmeta.n_groups_solute == 1463
    @test tmeta.n_groups_solvent == 14
    @test tmeta.unitcell ≈ [84.42188262939453, 84.42188262939453, 84.42188262939453]
end