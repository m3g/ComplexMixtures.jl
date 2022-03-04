
const bars = "-------------------------------------------------------------------------------"

"""

Print some information about the run.


"""
function title(R::Result, solute::Selection, solvent::Selection, nspawn::Int)

    println(bars)
    println("Starting MDDF calculation in parallel:")
    println("  $(R.nframes_read) frames will be considered")
    println("  Number of calculation threads: ", nspawn)
    println("  Solute: $(solute.natoms) atoms belonging to $(solute.nmols) molecules.")
    println("  Solvent: $(solvent.natoms) atom belonging to $(solvent.nmols) molecules.")

end

function title(R::Result, solute::Selection, solvent::Selection)

    println(bars)
    println("Starting MDDF calculation:")
    println("  $(R.nframes_read) frames will be considered")
    println("  Solute: $(solute.natoms) atoms belonging to $(solute.nmols) molecules.")
    println("  Solvent: $(solvent.natoms) atom belonging to $(solvent.nmols) molecules.")

end
