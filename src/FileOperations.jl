module FileOperations

#
# Function that determines the basename of a file, removing the path and the extension
#
clearname(filename::String) = remove_extension(basename(filename))

#
# Function that removes the extension of a file name
#
remove_extension(file::String) = file[1:findlast(==('.'), file)-1]

#
# Function that return only the extension of the file
#
file_extension(file::String) = file[findlast(==('.'), file)+1:end]

#
# Function that determines if a character is empty
#
empty_char(c::Char) = in(c,(Char(9),Char(32)))

#
# Function that checks if a line is a comment line or if it is empty
#
function commentary(s::String) 
    i = findfirst(c -> !(empty_char(c)), s) 
    return isnothing(i) || s[i] == '#'
end

end # module