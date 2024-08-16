"""
    write(
        R::Result, filename::String;
        solute_group_names::Vector{String} = R.solute.group_names,
        solvent_group_names::Vector{String} = R.solvent.group_names,
    )

Function to write the final results to output files as simple tables that are human-readable and easy to analyze with other software

If the solute and solvent group names are defined in `R`, the `solute_group_names` and `solvent_group_names` arguments are not necessary.
If they are not defined, the user can pass the names of the groups as strings in the `solute_group_names` and `solvent_group_names` arguments.

"""
function Base.write(
    R::Result, filename::String;
    solute_group_names::Union{Nothing,Vector{String}} = nothing,
    solvent_group_names::Union{Nothing,Vector{String}} = nothing,
)

    filename = expanduser(filename)

    # Set names of groups of solute and solvent (typically the atom names, but can
    # be any other name that the user wants to use, if custom groups were defined).
    solute_group_names = set_group_names(R, solute_group_names, :solute)
    solvent_group_names = set_group_names(R, solvent_group_names, :solvent)

    #
    # Total MDDF 
    #

    # Conversion factor for volumes (as KB integrals), from A^2 to cm^3/mol
    mole = 6.022140857e23
    convert = mole / 0.e24

    open(filename, "w") do output
        println(output, @sprintf("#"))
        println(output, "# Output of ComplexMixtures - MDDF - Version: $(VERSION)")
        println(output, @sprintf("#"))
        println(output, @sprintf("# Trajectory files and weights:"))
        for i in eachindex(R.files)
            println(output, "#  $(normpath(R.files[i].filename)) - w = $(R.weights[i])")
        end
        println(output, @sprintf("#"))
        println(output, @sprintf("# Density of solvent in simulation box (sites/A^2): %15.8f", R.density.solvent))
        println(output, @sprintf("# Density of solvent in bulk (estimated) (sites/A^2): %15.8f", R.density.solvent_bulk))
        println(output, @sprintf("# Molar volume of solvent in simulation (cc/mol): %14.8f", convert / R.density.solvent))
        println(output, @sprintf("# Molar volume of solvent in bulk (estimated) (cc/mol): %14.8f", convert / R.density.solvent_bulk))
        println(output, @sprintf("#"))
        println(output, @sprintf("# Number of atoms solute: %i8", natoms(R.solute)))
        println(output, @sprintf("# Number of atoms of the solvent: %i8", natoms(R.solvent)))
        println(output, @sprintf("#"))
        if R.files[1].options.usecutoff
            ibulk = setbin(R.files[1].options.dbulk, R.files[1].options.binstep)
            bulkerror = mean(R.mddf[ibulk:R.nbins])
            sdbulkerror = std(R.mddf[ibulk:R.nbins])
            println(output, "#")
            println(output, "# Using cutoff distance: $(R.cutoff)")
            println(output, @sprintf("# Average and standard deviation of bulk-gmd: %11.5f +/- %12.5f", bulkerror, sdbulkerror))
        end
        println(output, 
    """
    #
    # COLUMNS CORRESPOND TO:
    #       0  Minimum distance to solute (dmin)
    #       1  MDDF (md count normalized by md count of random-solute distribution)
    #       2  Kirwood-Buff integral (cc/mol) computed [(1/bulkdensity)*(col(6)-col(7))].
    #       3  Minimum distance site count for each dmin.
    #       4  Minimum distance site count for each dmin for random solute distribution.
    #       5  Cumulative number of molecules within dmin in the simulation (coordination number)
    #       6  Cumulative number of molecules within dmin for random solute distribution.
    #       7  Volume of the shell of distance dmin and width binstep.
    #
    #   0-DISTANCE         2-GMD      3-KB INT    4-MD COUNT  5-COUNT RAND      6-SUM MD    7-SUM RAND   8-SHELL VOL""")
        for i in eachindex(R.d)
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

    println("""

    OUTPUT FILES:

    Wrote total MDDF to output file: $filename
    """)

    # Solute and solvent contribution files
    solute_file = write_group_contributions(R, filename, :solute, solute_group_names) 
    solvent_file = write_group_contributions(R, filename, :solvent, solvent_group_names)

    return filename, solute_file, solvent_file
end

# Format numbers depending on their size
format(x) = abs(x) < 998 ? @sprintf("%12.7f", x) : @sprintf("%12.5e", x)

#
# Function that set the group names, given the optional input parameters of the write function
#
function set_group_names(R::Result, group_names::Union{Nothing,Vector{String}}, type::Symbol)
    atsel = getfield(R, type)
    atsel_str = type == :solute ? "Solute" : "Solvent"
    if isnothing(group_names)
        if !isempty(atsel.group_names)
            group_names = atsel.group_names
        else
            throw(ArgumentError("$atsel_str group names are not defined in R or as argument"))
        end
    else
        ngroups = isempty(atsel.group_atom_indices) ? natoms(atsel) : length(atsel.group_atom_indices)
        if ngroups != length(group_names)
            throw(ArgumentError("The number of $atsel_str group names is different from the number of groups."))
        end
    end
    return group_names
end

#
# Writting gmd per atom contributions for the solvent
#
function write_group_contributions(
    R::Result, 
    filename::String, 
    type::Symbol,
    group_names::Vector{String},
)
    if type == :solute
        atsel_str = "Solute"
        atsel = getfield(R, :solute)
        atsel_group_count = R.solute_group_count
    elseif type == :solvent
        atsel_str = "Solvent"
        atsel = getfield(R, :solvent)
        atsel_group_count = R.solvent_group_count
    else
        throw(ArgumentError("type must be either :solute or :solvent"))
    end

    # Names of output files containing atomic contibutions
    iext = findlast(==('.'), filename)
    if isnothing(iext)
        iext = length(filename) + 1
        extension = ".dat"
    else
        extension = filename[iext:end]
    end
    group_contributions_output = normpath(filename[1:iext-1] *
        "-GROUP_CONTRIBUTIONS_$(uppercase(atsel_str))" *
        extension
    )
    open(group_contributions_output, "w") do output
        println(output,
            chomp("""
            # $(atsel_str) atomic contributions to total MDDF.
            #
            # Trajectory files: 
            # $(join((file.filename for file in R.files), ",\n# "))
            #
            # Groups: 
            #
            """))
        for (i, name) in enumerate(group_names)
            println(output, @sprintf("# %8i %5s", i, name))
        end
        println(output, "#")
        # The first line of the table
        string = "#   DISTANCE    MDDF TOTAL"
        for i in eachindex(group_names)
            string = string * @sprintf("  %11i", i)
        end
        println(output, string)
        # The list of contributions
        for i in 1:R.nbins
            string = format(R.d[i])
            string = string * "  " * format(R.mddf[i]) # total mddf 
            if R.md_count_random[i] == 0
                string = string * "  " * repeat(format(0.0), length(group_names))
            else
                for group_count in atsel_group_count
                    contrib = group_count[i]
                    contrib = contrib / R.md_count_random[i]
                    string = string * "  " * format(contrib)
                end
            end
            println(output, string)
        end
    end # file writting

    println("Wrote solute group MDDF contributions to file: $group_contributions_output")

    return group_contributions_output
end

@testitem "write" begin
    using DelimitedFiles
    using ComplexMixtures
    using PDBTools
    using ComplexMixtures.Testing: data_dir
    atoms = readPDB("$data_dir/NAMD/structure.pdb")
    # Using or not bulk-range
    options1 = Options(stride = 5, seed = 321, StableRNG = true, nthreads = 1, silent = true, bulk_range=(8.0, 10.0))
    options2 = Options(stride = 5, seed = 321, StableRNG = true, nthreads = 1, silent = true, dbulk=8.0, usecutoff=false)
    protein = AtomSelection(select(atoms, "protein"), nmols = 1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol = 14)
    traj = Trajectory("$data_dir/NAMD/trajectory.dcd", protein, tmao, chemfiles = true)
    for options in (options1, options2)
        r = mddf(traj, options)
        tmpfile = tempname()*".dat"
        out1, out2, out3 = write(r, tmpfile)
        r_read = readdlm(out1, comments=true, comment_char='#')
        # Main output file
        @test r.d ≈ r_read[:,1]
        @test r.mddf ≈ r_read[:,2]
        @test r.kb ≈ r_read[:,3] rtol = 1e-5
        @test r.md_count ≈ r_read[:,4] 
        @test r.md_count_random ≈ r_read[:,5] rtol = 1e-5
        @test r.coordination_number ≈ r_read[:,6] rtol = 1e-5
        @test r.coordination_number_random ≈ r_read[:,7] rtol = 1e-5
        @test r.volume.shell ≈ r_read[:,8] rtol = 1e-5
        # Solute contributions
        r_read = readdlm(out2, comments=true, comment_char='#')
        @test r.d ≈ r_read[:,1]
        @test r.mddf ≈ r_read[:,2] 
        for i in eachindex(r.solute.group_names)
            @test contributions(r, SoluteGroup([i])) ≈ r_read[:,i+2] rtol = 1e-5
        end
        # Solvent contributions
        r_read = readdlm(out3, comments=true, comment_char='#')
        @test r.d ≈ r_read[:,1]
        @test r.mddf ≈ r_read[:,2] 
        for (i, name) in enumerate(atom_group_names(r.solvent))
            @test contributions(r, SolventGroup([name])) ≈ r_read[:,i+2] rtol = 1e-5
        end
    end
end
