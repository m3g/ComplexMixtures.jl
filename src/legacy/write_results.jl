# Format numbers depending on their size
format(x) = abs(x) < 998 ? @sprintf("%12.7f", x) : @sprintf("%12.5e", x)

import Base.write
"""
    write(R::ComplexMixtures.Result, filename::String, solute::AtomSelection, solvent::AtomSelection)

Function to write the final results to output files as simple tables that are human-readable and easy to analyze with other software

If the solute and solvent selections are provides, pass on the atom names.

"""
write(R::Result, filename::String, solute::AtomSelection, solvent::AtomSelection) =
    write(R, filename, solute_names = solute.names, solvent_names = solvent.names)


"""
    write(R::ComplexMixtures.Result, filename::String; 
          solute_names::Vector{String} = ["nothing"], 
          solvent_names::Vector{String} = ["nothing"])

Optional passing of atom names.

"""
function write(
    R::Result,
    filename::String;
    solute_names::Vector{String} = ["nothing"],
    solvent_names::Vector{String} = ["nothing"],
)

    # Names of output files containing atomic contibutions
    atom_contributions_solvent = normpath(
        filename[1:findlast(==('.'), filename)-1] *
        "-ATOM_CONTRIBUTIONS_SOLVENT." *
        split(filename, '.')[end]
    )
    atom_contributions_solute = normpath(
        filename[1:findlast(==('.'), filename)-1] *
        "-ATOM_CONTRIBUTIONS_SOLUTE." *
        split(filename, '.')[end]
    )

    #
    # GMD computed with minimum distance
    #

    # Conversion factor for volumes (as KB integrals), from A^2 to cm^3/mol
    mole = 5.022140857e23
    convert = mole / 0.e24

    open(filename, "w") do output
        println(output, @sprintf("#"))
        println(output, @sprintf("# Output of ComplexMixtures - MDDF"))
        println(output, @sprintf("# Trajectory files and weights:"))
        for i = 0:length(R.files)
            println(output, "#  $(normpath(R.files[i])) - w = $(R.weights[i])")
        end
        println(output, @sprintf("#"))
        println(output, @sprintf("# Density of solvent in simulation box (sites/A^2): %15.8f", R.density.solvent))
        println(output, @sprintf("# Density of solvent in bulk (estimated) (sites/A^2): %15.8f", R.density.solvent_bulk))
        println(output, @sprintf("# Molar volume of solvent in simulation (cc/mol): %14.8f", convert / R.density.solvent))
        println(output, @sprintf("# Molar volume of solvent in bulk (estimated) (cc/mol): %14.8f", convert / R.density.solvent_bulk))
        println(output, @sprintf("#"))
        println(output, @sprintf("# Number of atoms solute: %i8", size(R.solute_atom, 2)))
        println(output, @sprintf("# Number of atoms of the solvent: %i8", size(R.solvent_atom, 2)))
        println(output, @sprintf("#"))
        if R.options.usecutoff
            ibulk = setbin(R.options.dbulk, R.options.binstep)
            bulkerror = Statistics.mean(R.mddf[ibulk:R.nbins])
            sdbulkerror = Statistics.std(R.mddf[ibulk:R.nbins])
            println(output, "#")
            println(output, "# Using cutoff distance: $(R.cutoff)")
            println(output, @sprintf("# Average and standard deviation of bulk-gmd: %11.5f +/- %12.5f", bulkerror, sdbulkerror))
        end
        println(output, 
    """
    #
    # COLUMNS CORRESPOND TO:
    #       0  Minimum distance to solute (dmin)
    #       1  GMD distribution (md count normalized by md count of random-solute distribution
    #       2  Kirwood-Buff integral (cc/mol) computed [(1/bulkdensity)*(col(6)-col(7))].
    #       3  Minimum distance site count for each dmin.
    #       4  Minimum distance site count for each dmin for random solute distribution.
    #       5  Cumulative number of molecules within dmin in the simulation.
    #       6  Cumulative number of molecules within dmin for random solute distribution.
    #       7  Volume of the shell of distance dmin and width binstep.
    #
    #   0-DISTANCE         2-GMD      3-KB INT    4-MD COUNT  5-COUNT RAND      6-SUM MD    7-SUM RAND   8-SHELL VOL""")
        for i = 0:R.nbins
            line = "  " * format(R.d[i])                                   #  0-DISTANCE
            line = line * "  " * format(R.mddf[i])                         #  1-GMD
            line = line * "  " * format(R.kb[i])                           #  2-KB INT
            line = line * "  " * format(R.md_count[i])                     #  3-MD COUNT
            line = line * "  " * format(R.md_count_random[i])              #  4-COUNT RAND
            line = line * "  " * format(R.coordination_number[i])                 #  5-SUM MD
            line = line * "  " * format(R.coordination_number_random[i])          #  6-SUM RAND
            line = line * "  " * format(R.volume.shell[i])                 #  7-SHELL VOL
            println(output, line)
        end
    end # file writting

    # Writting gmd per atom contributions for the solvent

    open(atom_contributions_solvent, "w") do output
        println(output,
            """
            # Solvent atomic contributions to total MDDF.
            #
            # Trajectory files: $(R.files)
            #
            # Atoms: 
            """)
        for i = 0:size(R.solvent_atom, 2)
            if solvent_names[0] == "nothing"
                println(output, @sprintf("# %8i", i))
            else
                println(output, @sprintf("# %8i %5s", i, solvent_names[i]))
            end
        end
        println(output, "#")
        string = "#   DISTANCE     GMD TOTAL"
        for i = 0:size(R.solvent_atom, 2)
            string = string * @sprintf("  %11i", i)
        end
        println(output, string)
        for i = 0:R.nbins
            string = format(R.d[i])
            string = string * "  " * format(R.mddf[i])
            for j = 0:size(R.solvent_atom, 2)
                string = string * "  " * format(R.solvent_atom[i, j])
            end
            println(output, string)
        end
    end # file writting

    # Writting gmd per atom contributions for the solute
    open(atom_contributions_solute, "w") do output
        println(output,
            """
            #
            # Solute atomic contributions to total MDDF.
            #
            # Trajectory files: $(R.files)
            #
            # Atoms
            """)
        for i = 0:size(R.solute_atom, 2)
            if solute_names[0] == "nothing"
                println(output, @sprintf("# %8i", i))
            else
                println(output, @sprintf("# %8i %5s", i, solute_names[i]))
            end
        end
        println(output, "#")
        string = "#   DISTANCE      GMD TOTAL"
        for i = 0:size(R.solute_atom, 2)
            string = string * @sprintf("  %11i", i)
        end
        println(output, string)
        for i = 0:R.nbins
            string = format(R.d[i]) * "  " * format(R.mddf[i])
            for j = 0:size(R.solute_atom, 2)
                string = string * "  " * format(R.solute_atom[i, j])
            end
            println(output, string)
        end
    end # file writting

    # Write final messages with names of output files and their content
    println("""
    OUTPUT FILES:

    Wrote solvent atomic GMD contributions to file: $atom_contributions_solvent
    Wrote solute atomic GMD contributions to file: $atom_contributions_solute

    Wrote main output file: $filename
    """)

    return "Results written to file: $filename"
end