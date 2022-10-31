"""
    Trajectory(filename::String, solute::Selection, solvent::Selection; format::String = "", chemfiles = false)

Trajectory constructor data type. 

Defaults to reading with the Chemfiles infrastructure, except for DCD and PDB trajectory
files, if the "PDBTraj" option is provided.

See memory issue (https://github.com/chemfiles/Chemfiles.jl/issues/44)

"""
abstract type Trajectory end

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
