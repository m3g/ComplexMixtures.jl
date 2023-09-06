"""
    Trajectory(filename::String, solute::Selection, solvent::Selection; format::String = "", chemfiles = false)

Trajectory constructor data type. 

Defaults to reading with the Chemfiles infrastructure, except for DCD and PDB trajectory
files, if the "PDBTraj" option is provided.

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

# Include specific predefined trajectory types
include("./trajectory_formats/ChemFiles.jl")
include("./trajectory_formats/NamdDCD.jl")
include("./trajectory_formats/PDBTraj.jl")

function Trajectory(
    filename::String,
    solute::Selection,
    solvent::Selection;
    format::String = "",
    chemfiles = false,
)
    if !chemfiles && (format == "dcd" || FileOperations.file_extension(filename) == "dcd")
        trajectory = NamdDCD(filename, solute, solvent)
    elseif !chemfiles && format == "PDBTraj"
        trajectory = PDBTraj(filename, solute, solvent)
    else
        trajectory = ChemFile(filename, solute, solvent, format = format)
    end
    return trajectory
end

# If only one selection is provided, assume that the solute and the solvent are the same
Trajectory(filename::String, solvent::Selection; format::String = "", chemfiles = false) =
    Trajectory(filename, solvent, solvent, format = format, chemfiles = chemfiles)

#=
    convert_unitcell(unitcell::Union{SVector{3}, SMatrix{3,3}})

Function to return the unit cell as a vector or matrix, depending on if 
the cell is diagonal or not, up to a relative precision of 1e-10 by default.

=#
function convert_unitcell(unitcell::AbstractMatrix; tol = 1e-10)
    size(unitcell) == (3,3) || error("Unit cell must be a 3x3 matrix.")
    s = minimum(diag(unitcell))
    is_diag = all(unitcell[i,j] < tol*s for i in 1:3, j in 1:3 if i != j)
    return is_diag ?  SVector(diag(unitcell)) : SMatrix(unitcell)
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
    protein = Selection(select(atoms, "protein"), nmols = 1)
    tmao = Selection(select(atoms, "resname TMAO"), natomspermol = 14)

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
        format = "PDBTraj",
    )
    @test traj.solute.natoms == 1463
    @test traj.solvent.natoms == 2534
    @test ComplexMixtures.convert_unitcell(ComplexMixtures.getunitcell(traj)) ≈ SVector(84.47962951660156, 84.47962951660156, 84.47962951660156)

    # Chemfiles with NAMD
    traj = Trajectory(
        "$(Testing.data_dir)/NAMD/trajectory.dcd",
        protein,
        tmao,
        chemfiles = true,
    )
    @test traj.nframes == 20
    @test ComplexMixtures.convert_unitcell(ComplexMixtures.getunitcell(traj)) ≈ SVector(84.42188262939453, 84.42188262939453, 84.42188262939453)
    @test traj.solute.natoms == 1463
    @test traj.solvent.natoms == 2534

    # Chemfiles with Gromacs
    atoms = readPDB("$(Testing.data_dir)/Gromacs/system.pdb")
    protein = Selection(select(atoms, "protein"), nmols = 1)
    emi = Selection(select(atoms, "resname EMI"), natomspermol = 20)
    traj = Trajectory("$(Testing.data_dir)/Gromacs/trajectory.xtc", protein, emi)
    @test traj.nframes == 26
    @test ComplexMixtures.convert_unitcell(ComplexMixtures.getunitcell(traj)) ≈ SVector(95.11481285095215, 95.11481285095215, 95.13440132141113)
    @test traj.solute.natoms == 1231
    @test traj.solvent.natoms == 5080

end
