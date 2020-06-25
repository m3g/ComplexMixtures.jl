#
# This is just a wrapper to the searchsortedfirst function in base
# which returns 0 if the elements is not found in the list instead
# of the position of the putative insertion
#
function my_searchsortedfirst(x,i)
  index = searchsortedfirst(x,i)
  if index > lastindex(x) || x[index] != i
    return 0
  else
    return index
  end
end
