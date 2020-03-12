# 
# Select atoms using vmd selection syntax, with vmd in background
#

# Function to return the selection from a input file (topology, coordinates, etc), 
# by calling VMD in the background.

function VMDselect( inputfile :: String, selection :: String; vmd="vmd" )

  local index_list :: String 
  local readnext :: Bool = false

  vmd_input = Base.open("./VMDINPUT_TMP.VMD","w")
  Base.write(vmd_input,"mol new \"$inputfile\" \n")
  Base.write(vmd_input,"set sel [ atomselect top \"$selection\" ] \n")
  Base.write(vmd_input,"puts \"INDEXLIST\" \n")
  Base.write(vmd_input,"set indexes [ \$sel get index ] \n")
  Base.write(vmd_input,"puts \"ENDINDEXLIST\" \n")
  Base.write(vmd_input,"exit \n")
  Base.close(vmd_input)

  vmd_output = read(`$vmd -dispdev text -e ./VMDINPUT_TMP.VMD`, String)

  for line in split(vmd_output,"\n")
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
  selection_indexes = Vector{Int64}(undef,nsel) 
  for i in 1:nsel
    selection_indexes[i] = parse(Int64,index_split[i]) + 1
  end

  run(`\rm -f ./VMDINPUT_TMP.VMD`)

  return selection_indexes

end














