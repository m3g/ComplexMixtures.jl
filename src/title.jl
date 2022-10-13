
const bars = "-------------------------------------------------------------------------------"

atoms_str(n) = "$n $(n == 1 ? "atom" : "atoms")"
mol_str(n) = "$n $(n == 1 ? "molecule" : "molecules")"

"""

Print some information about the run.


"""
function title(R::Result, solute::Selection, solvent::Selection, nspawn::Int)
    println(bars)
    println("Starting MDDF calculation in parallel:")
    println("  $(R.nframes_read) frames will be considered")
    println("  Number of calculation threads: ", nspawn)
    println("  Solute: $(atoms_str(solute.natoms)) belonging to $(mol_str(solute.nmols)).")
    println("  Solvent: $(atoms_str(solvent.natoms)) belonging to $(mol_str(solvent.nmols)).")

end

function title(R::Result, solute::Selection, solvent::Selection)
    println(bars)
    println("Starting MDDF calculation:")
    println("  $(R.nframes_read) frames will be considered")
    println("  Solute: $(atoms_str(solute.natoms)) belonging to $(mol_str(solute.nmols)).")
    println("  Solvent: $(atoms_str(solvent.natoms)) belonging to $(mol_str(solvent.nmols))")

end
