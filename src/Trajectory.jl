# Default reading with the Chemfiles infrastructure

Trajectory( filename :: String, 
            solute :: Selection, solvent :: Selection;
            format :: String = "") =
  ChemFile(filename,solute,solvent,format=format)

# If only one selection is provided, assume that the solute and the 
# solvent are the same

Trajectory( filename :: String, 
            solvent :: Selection; 
            format :: String = "") =
  ChemFile(filename,solvent,solvent,format=format)

