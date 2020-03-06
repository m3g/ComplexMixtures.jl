module FileOperations

  #
  # Function that determines the basename of a file,
  # removing the path and the extension
  #

  clearname(filename :: String) = remove_extension(basename(filename))
  
  #
  # Function that removes the extension of a file name
  #
  
  function remove_extension(filename :: String)
    i = length(filename)
    idot = i+1
    while i > 0
      if filename[i] == '.' 
        idot = i
        break
      end
      i = i - 1
    end
    return filename[1:idot-1]
  end

  #
  # Function that return only the extension of the file
  #
  
  function file_extension(filename :: String)
    i = length(filename)
    idot = 1
    while i > 0 
      if filename[i] == '.'
        idot = i
        break
      end
      i = i - 1
    end
    return filename[idot+1:length(filename)]
  end

  #
  # Function that determines if a character is empty
  #
  
  function empty_char(char :: Char)
    if char == Char(9)  || char == Char(32)
      return true
    else
      return false
    end
  end
  
  #
  # Function that checks if a line is a comment line or if it is empty
  #

  function commentary(string :: String)
    i = 1
    while empty_char(string[i]) && i < length(string)
      i = i + 1
    end
    if string[i] == '#' || i == length(string) 
      return true
    else
      return false
    end
  end

end
