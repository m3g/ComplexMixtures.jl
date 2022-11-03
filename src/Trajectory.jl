"""
    Trajectory(filename::String, solute::Selection, solvent::Selection; format::String = "", chemfiles = false)

Trajectory constructor data type. 

Defaults to reading with the Chemfiles infrastructure, except for DCD and PDB trajectory
files, if the "PDBTraj" option is provided.

See memory issue (https://github.com/chemfiles/Chemfiles.jl/issues/44)
"""
abstract type Trajectory end

include("./trajectory_formats/ChemFiles.jl")
include("./trajectory_formats/NamdDCD.jl")
include("./trajectory_formats/PDBTraj.jl")

function Trajectory(filename::String, solute::Selection, solvent::Selection; format::String = "", chemfiles = false)
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

@testitem "Trajectory" begin
    using ComplexMixtures
    using ComplexMixtures.Testing
    using PDBTools
    using StaticArrays

    atoms = readPDB(Testing.pdbfile)  
    protein = Selection(select(atoms,"protein"),nmols=1)
    tmao = Selection(select(atoms,"resname TMAO"),natomspermol=14)

    # NAMD DCD file
    traj = Trajectory("$(Testing.data_dir)/NAMD/trajectory.dcd",protein,tmao)
    @test traj.nframes == 20
    @test traj.lastatom == 4012
    @test traj.sides_in_dcd == true

    # PDB file
    traj = Trajectory("$(Testing.data_dir)/PDB/trajectory.pdb",protein,tmao,format="PDBTraj")
    @test traj.natoms == 62026
    @test traj.solute.natoms == 1463
    @test traj.solvent.natoms == 2534
    @test traj.sides[1] ≈ SVector(84.42188262939453, 84.42188262939453, 84.42188262939453)

    # Chemfiles with NAMD
    traj = Trajectory("$(Testing.data_dir)/NAMD/trajectory.dcd",protein,tmao,chemfiles=true)
    @test traj.natoms == 62026
    @test traj.nframes == 20
    @test traj.sides[1] ≈ SVector(84.42188262939453, 84.42188262939453, 84.42188262939453)
    @test traj.solute.natoms == 1463
    @test traj.solvent.natoms == 2534
    
    # Chemfiles with Gromacs
    atoms = readPDB("$(Testing.data_dir)/Gromacs/system.pdb")  
    protein = Selection(select(atoms,"protein"),nmols=1)
    emi = Selection(select(atoms,"resname EMI"),natomspermol=20)
    traj = Trajectory("$(Testing.data_dir)/Gromacs/trajectory.xtc",protein,emi)
    @test traj.natoms == 84864
    @test traj.nframes == 26
    @test traj.sides[1] ≈ SVector(95.11481285095215, 95.11481285095215, 95.13440132141113)
    @test traj.solute.natoms == 1231
    @test traj.solvent.natoms == 5080

end
