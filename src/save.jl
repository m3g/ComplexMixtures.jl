using StructTypes
using JSON3

StructTypes.StructType(::Type{SolSummary}) = StructTypes.Struct()
StructTypes.StructType(::Type{MutableResult}) = StructTypes.Struct()
StructTypes.StructType(::Type{Result}) = StructTypes.Struct()
StructTypes.StructType(::Type{Density}) = StructTypes.Struct()
StructTypes.StructType(::Type{Volume}) = StructTypes.Struct()
StructTypes.StructType(::Type{Options}) = StructTypes.Struct()

#
# Function to write the result data structure to a json file
#

function save(R::Result, filename::String)
  f = open(filename,"w")
  JSON3.write(f,R)
  close(f)
end

