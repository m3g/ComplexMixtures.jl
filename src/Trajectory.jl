#
# Default reading with the Chemfiles infrastructure, except for DCD files,
# because of the memory issue (https://github.com/chemfiles/Chemfiles.jl/issues/44)
#

function Trajectory( filename :: String, 
                     solute :: Selection, solvent :: Selection;
                     format :: String = "")
  if format == "dcd" || FileOperations.file_extension(filename) == "dcd"
    NamdDCD(filename,solute,solvent)
  else
    ChemFile(filename,solute,solvent,format=format)
  end
end

# If only one selection is provided, assume that the solute and the 
# solvent are the same

function Trajectory( filename :: String, 
                    solvent :: Selection; 
                    format :: String = "")
  if format == "dcd" || FileOperations.file_extension(filename) == "dcd"
    NamdDCD(filename,solvent,solvent)
  else
    ChemFile(filename,solvent,solvent,format=format)
  end
end

