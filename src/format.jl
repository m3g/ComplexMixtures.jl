# Function to format propertly output numbers in tables

function format(x)

  if abs(x) < 999.
    format = @sprintf("%12.7f",x)
  else
    format = @sprintf("%12.5f",x)
  end
  return format

end

