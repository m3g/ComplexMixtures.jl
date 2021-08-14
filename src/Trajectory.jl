"""

```
Trajectory
```

Trajectory data type. 
Default reading with the Chemfiles infrastructure, except for DCD and PDB trajectory
files, if the "PDBTraj" option is provided.
See memory issue (https://github.com/chemfiles/Chemfiles.jl/issues/44)

"""
abstract type Trajectory end

function Trajectory(filename::String,
                    solute::Selection, solvent::Selection;
                    format::String = "",
                    chemfiles=false)
  if !chemfiles 
    if format == "dcd" || FileOperations.file_extension(filename) == "dcd"
      NamdDCD(filename,solute,solvent)
    elseif format == "PDBTraj"
      PDBTraj(filename,solute,solvent)
    else
      ChemFile(filename,solute,solvent,format=format)
    end
  else
    ChemFile(filename,solute,solvent,format=format)
  end
end

# If only one selection is provided, assume that the solute and the 
# solvent are the same

function Trajectory(filename::String,
                    solvent::Selection;
                    format::String = "",
                    chemfiles=false)
  Trajectory(filename,solvent,solvent,format=format,chemfiles=chemfiles)
end

