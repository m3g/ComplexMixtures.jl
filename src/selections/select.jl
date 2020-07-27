#
# This is a simple extension to the selection function of Chemfiles, which makes the 
# selection one-indexed, returns the indexes in Int64 (not UInt64) format, and adds 
# some simple selection operators, such as "protein"
#
# Usually the "filename" here is a PDB or topology file within the formats accepted by Chemfiles
#

import Chemfiles

function select(filename :: String, selection :: String ; format="", onebased=true)

  # Parse some additional useful types of selections
  selection = select_macros(selection)

  # Convert the user input to zero based, that is understood by Chemfiles
  if onebased
    selection = index_to_zero_based(selection)
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

# If there is a selection by indexes, automatically replace the indexes by the zero-based
# indexing to pass it to Chemfiles

function index_to_zero_based(selection)
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
  end
  selection = ""
  for item in s
    selection = selection*"$item "
  end
  return rstrip(selection)
end


