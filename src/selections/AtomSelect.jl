#
# This is a simple extension to the selection function of Chemfiles, which makes the 
# selection one-indexed, returns the indexes in Int64 (not UInt64) format, and adds 
# some simple selection operators, such as "protein"
#
# Usually the "filename" here is a PDB or topology file within the formats accepted by Chemfiles
#

module AtomSelect

  import Chemfiles
  
  #
  # The final easy-to-use select function
  #
  function select(filename :: String, selection :: String ; format="", onebased=true)
  
    # Parse some additional useful types of selections
    selection = select_macros(selection)
  
    # Convert the user input to zero based, that is understood by Chemfiles
    if onebased
      selection = indexes_to_zero_based(selection)
    end
  
    # Open the trajectory, define the selection, get frame and indexes, using Chemfiles
    t = Chemfiles.Trajectory(filename,'r',format)
    sel = Chemfiles.Selection(selection)
    frame = Chemfiles.read(t)
    indexes = Chemfiles.evaluate(sel,frame)
  
    # Convert indexes back to one-based and to Int type
    if onebased
      for i in 1:length(indexes)
        indexes[i] = (indexes[i]%Int)+1
      end
    end
    Chemfiles.close(t)
  
    return indexes
  end

  # 
  # Convert indexes to zero-based indexing to pass to Chemfiles
  #
  function indexes_to_zero_based(selection :: String)
    # 'Index' keyword
    if occursin("index",selection)
      s = split(selection)
      ikey = 1
      while ikey <= length(s)  
        if s[ikey] == "index"
          ikey += 1
          if tryparse(Int,s[ikey]) == nothing
            ikey += 1
          end
          index = tryparse(Int,s[ikey])
          while ikey <= length(s) && index != nothing 
            index = index - 1
            s[ikey] = "$index"
            ikey += 1
            if ikey <= length(s)
              index = tryparse(Int,s[ikey])
            end
          end
        end
        ikey += 1
      end
      selection = ""
      for item in s
        selection = selection*"$item "
      end
      selection = rstrip(selection)
    end
    return selection
  end
  
  #
  # Some macro-selection missing in chemfiles that are useful for simulations
  #
  function select_macros(selection :: String)
     macro_list = [ 
                    [ "protein", "(resname ALA ARG ASN ASP ASX CYS GLU GLN GLX GLY HIS HSD ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL)" ],
                    [ "water", "(resname TIP3P TIP3 TIP4P TIP4 HOH WAT OH2 SPCE SPC TIP5P TIP7P)" ]
                  ]
     # Select all that is protein
     for item in macro_list 
       if occursin("$(item[1])",selection)
         selection = replace(selection, "$(item[1])" => "$(item[2])")
       end
     end
  
     return selection
  end

  # 
  # Convert indexes inside function calls to zero-based as well
  #
  function functions_to_zero_based(selection :: String)
    #
    # List of functions that use numeric indexes and number of arguments
    #
    funcs = [ 
              [ "is_bonded", 2 ],
              [ "IS_ANGLE", 3 ], 
              [ "IS_DIHEDRAL", 4 ],
              [ "is_improper", 4 ],
              [ "distance", 2 ],
              [ "angle", 3 ],
              [ "dihedral", 4 ],
              [ "out_of_plane", 4 ]
            ]
    # First, disambiguate is_angle and is_diehdral 
    selection = replace(selection, "is_angle" => "IS_ANGLE")
    selection = replace(selection, "is_dihedral" => "IS_DIHEDRAL")
    # Cycle over the possibe functions that receive indexes as arguments 
    for f in funcs
      ranges = findall(f[1],selection)
      l = length(f[1])
      for range in ranges
        # find closing parenthesis
        open = 1
        pos = range[1]+l
        while open > 0
          pos += 1
          if selection[pos:pos] == ")"
            open -= 1
          elseif selection[pos:pos] == "("
            open += 1
          end
          if open > 1 && selection[pos:pos] == ","
            new_selection = selection[1:pos-1]*"@"*selection[pos+1:length(selection)]
            selection = new_selection
          end
        end
        # Split the arguments
        args = split(selection[range[1]+l+1:pos-1],',')
        # Try to parse each of the arguments into an integer, and if
        # so subtract 1 
        for iarg in 1:f[2]
          index = tryparse(Int,args[iarg])
          if index != nothing
            index -= 1
            args[iarg] = "$index"
          end
        end
        # Replace the arguments of the function with the new indexes
        new_selection = selection[1:range[l]+1]
        for iarg in 1:f[2]-1
          new_selection = new_selection*args[iarg]*","
        end
        new_selection = new_selection*args[f[2]]
        selection = new_selection*selection[pos:length(selection)]
        selection = replace(selection,"@" => ",")
      end
    end
    # Restore function names
    selection = replace(selection, "IS_ANGLE" => "is_angle")
    selection = replace(selection, "IS_DIHEDRAL" => "is_dihedral")
    return selection
  end

  export select
end

