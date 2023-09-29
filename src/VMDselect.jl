""" 
    VMDselect(inputfile::String, selection::String; vmd="vmd", srcload=nothing)

Select atoms using vmd selection syntax, with vmd in background

Returns the list of index (one-based) and atom names

Function to return the selection from a input file (topology, coordinates, etc), 
by calling VMD in the background.

The `srcload` argument can be used to load a list of scripts before loading the input file,
for example with macros to define custom selection keywords.

"""
function VMDselect(inputfile::String, selection::String; vmd = "vmd", srcload = nothing)

    if !isfile(inputfile)
        error("Could not find file: $inputfile")
    end

    vmdinput_file = tempname()
    vmd_input = Base.open(vmdinput_file, "w")
    if !isnothing(srcload)
        if srcload isa AbstractString
            srcload = [srcload]
        end
        for srcfile in srcload
            Base.write(vmd_input, "source \"$srcfile\" \n")
        end
    end
    Base.write(vmd_input, "mol new \"$inputfile\" \n")
    Base.write(vmd_input, "set sel [ atomselect top \"$selection\" ] \n")
    Base.write(vmd_input, "puts \"INDEXLIST\" \n")
    Base.write(vmd_input, "set indexes [ \$sel get index ] \n")
    Base.write(vmd_input, "puts \"ENDINDEXLIST\" \n")
    Base.write(vmd_input, "puts \"NAMELIST\" \n")
    Base.write(vmd_input, "set names [ \$sel get name ] \n")
    Base.write(vmd_input, "puts \"ENDNAMELIST\" \n")
    Base.write(vmd_input, "exit \n")
    Base.close(vmd_input)

    vmd_output = Base.read(`$vmd -dispdev text -e $vmdinput_file`, String)

    # Read indexes
    local index_list::String
    readnext = false
    for line in split(vmd_output, "\n")
        if readnext
            if line == "ENDINDEXLIST"
                error("ERROR: Selection '$selection' does not contain any atom")
            end
            index_list = line
            break
        end
        if line == "INDEXLIST"
            readnext = true
        end
    end
    index_split = split(index_list)
    nsel = length(index_split)
    selection_indexes = Vector{Int}(undef, nsel)
    for i = 1:nsel
        selection_indexes[i] = parse(Int, index_split[i]) + 1
    end

    # Read atom names
    local name_list::String
    readnext = false
    for line in split(vmd_output, "\n")
        if readnext
            if line == "ENDNAMELIST"
                error("ERROR: Selection '$selection' does not contain any atom")
            end
            name_list = line
            break
        end
        if line == "NAMELIST"
            readnext = true
        end
    end
    name_split = split(name_list)
    nsel = length(name_split)
    selection_names = Vector{String}(undef, nsel)
    for i = 1:nsel
        selection_names[i] = strip(name_split[i])
    end

    return selection_indexes, selection_names
end

@testitem "VMDSelect" begin
    using ComplexMixtures
    pdbfile = "$(@__DIR__)/../test/data/NAMD/structure.pdb"
    if !isnothing(Sys.which("vmd"))
        @test VMDselect(pdbfile, "protein and residue 1") == (
            [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23],
            ["N", "HN", "CA", "HA", "CB", "HB1", "HB2", "SG", "HG1", "C", "O"],
        )
    end
end
