#
# Function that checks if the distance is in bulk or not
# according to the use of cutoff or not
#
function inbulk(d :: Float64, R :: Result)
  if ! R.options.usecutoff
    d > R.dbulk
  else
    d > R.dbulk && d <= R.cutoff   
  end
end
