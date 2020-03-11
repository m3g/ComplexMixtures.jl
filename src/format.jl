# Function to format propertly output numbers in tables

function format(x)

  if abs(x) < 999.
    return @sprintf("%12.7f",x)
  else
    return @sprintf("%12.5f",x)
  end

end

