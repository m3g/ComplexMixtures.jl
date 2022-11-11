"""
    Trajectory(filename::String, solute::Selection, solvent::Selection; format::String = "", chemfiles = false)

Trajectory constructor data type. 

Defaults to reading with the Chemfiles infrastructure, except for DCD and PDB trajectory
files, if the "PDBTraj" option is provided.

See memory issue (https://github.com/chemfiles/Chemfiles.jl/issues/44)
"""
abstract type Trajectory end

# Include specific predefined trajectory types
include("./trajectory_formats/ChemFiles.jl")
include("./trajectory_formats/NamdDCD.jl")
include("./trajectory_formats/PDBTraj.jl")

function Trajectory(filename::String, solute::Selection, solvent::Selection; format::String="", chemfiles=false)
    if !chemfiles && (format == "dcd" || FileOperations.file_extension(filename) == "dcd")
        trajectory = NamdDCD(filename, solute, solvent)
    elseif !chemfiles && format == "PDBTraj"
        trajectory = PDBTraj(filename, solute, solvent)
    else
        trajectory = ChemFile(filename, solute, solvent, format=format)
    end
    return trajectory
end

# If only one selection is provided, assume that the solute and the solvent are the same
Trajectory(filename::String, solvent::Selection; format::String="", chemfiles=false) =
    Trajectory(filename, solvent, solvent, format=format, chemfiles=chemfiles)

#
# Function to get the appropriate representation of the unit cell, depending on its type
#
function setunitcell(uc::AbstractVecOrMat)
    unitcell = if uc isa AbstractVector # Orthorhombic cell
        SVector(uc)
    else # Orthorhombic in practice 
        if isdiag(uc)
            SVector{3}(uc[i,i] for i = 1:3)
        else # Triclinic cell
            SMatrix{3,3}(uc)
        end
    end
    return unitcell
end
# From the trajectory file
setunitcell(trajectory::Trajectory) = setunitcell(trajectory.sides[1])

@testitem "setunitcell" begin
    using ComplexMixtures
    using ComplexMixtures.Testing
    using PDBTools
    using StaticArrays

    # From the definition of the unitcell
    uc = SVector(1.0, 1.0, 1.0)
    @test ComplexMixtures.setunitcell(uc) == uc
    uc_mat = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    @test ComplexMixtures.setunitcell(uc_mat) == uc
    uc_mat = SMatrix{3,3}(1.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    @test ComplexMixtures.setunitcell(uc_mat) == uc_mat

    # From the trajectory
    atoms = readPDB(Testing.pdbfile)
    options = Options(stride=5,seed=321,StableRNG=true,nthreads=1,silent=true)
    protein = Selection(select(atoms, "protein"), nmols=1)
    tmao = Selection(select(atoms, "resname TMAO"), natomspermol=14)
    traj = Trajectory("$(Testing.data_dir)/NAMD/trajectory.dcd", protein, tmao)
    ComplexMixtures.nextframe!(traj)
    @test ComplexMixtures.setunitcell(traj) == SVector(84.42188262939453, 84.42188262939453, 84.42188262939453)
    ComplexMixtures.closetraj(traj)
end

@testitem "Trajectory" begin
    using ComplexMixtures
    using ComplexMixtures.Testing
    using PDBTools
    using StaticArrays

    atoms = readPDB(Testing.pdbfile)
    protein = Selection(select(atoms, "protein"), nmols=1)
    tmao = Selection(select(atoms, "resname TMAO"), natomspermol=14)

    # NAMD DCD file
    traj = Trajectory("$(Testing.data_dir)/NAMD/trajectory.dcd", protein, tmao)
    @test traj.nframes == 20
    @test traj.lastatom == 4012
    @test traj.sides_in_dcd == true

    # PDB file
    traj = Trajectory("$(Testing.data_dir)/PDB/trajectory.pdb", protein, tmao, format="PDBTraj")
    @test traj.natoms == 62026
    @test traj.solute.natoms == 1463
    @test traj.solvent.natoms == 2534
    @test traj.sides[1] ≈ SVector(84.42188262939453, 84.42188262939453, 84.42188262939453)

    # Chemfiles with NAMD
    traj = Trajectory("$(Testing.data_dir)/NAMD/trajectory.dcd", protein, tmao, chemfiles=true)
    @test traj.natoms == 62026
    @test traj.nframes == 20
    @test traj.sides[1] ≈ SVector(84.42188262939453, 84.42188262939453, 84.42188262939453)
    @test traj.solute.natoms == 1463
    @test traj.solvent.natoms == 2534

    # Chemfiles with Gromacs
    atoms = readPDB("$(Testing.data_dir)/Gromacs/system.pdb")
    protein = Selection(select(atoms, "protein"), nmols=1)
    emi = Selection(select(atoms, "resname EMI"), natomspermol=20)
    traj = Trajectory("$(Testing.data_dir)/Gromacs/trajectory.xtc", protein, emi)
    @test traj.natoms == 84864
    @test traj.nframes == 26
    @test traj.sides[1] ≈ SVector(95.11481285095215, 95.11481285095215, 95.13440132141113)
    @test traj.solute.natoms == 1231
    @test traj.solvent.natoms == 5080

end
