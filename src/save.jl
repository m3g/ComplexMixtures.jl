using StructTypes
using JSON3

StructTypes.StructType(::Type{Result}) = StructTypes.Struct()
StructTypes.StructType(::Type{Density}) = StructTypes.Struct()
StructTypes.StructType(::Type{Volume}) = StructTypes.Struct()
StructTypes.StructType(::Type{Options}) = StructTypes.Struct()

function write( R :: Result, filename :: String)
  f = open(filename,"w")
  JSON3.write(f,R)
  close(f)
end

function read(filename :: String)
  R = read_wg(filename)
  GC.gc()
  return R
end

function read_wg(filename :: String)
  f = open(filename,"r")
  R = JSON3.read(f,Result)
  # Need to convert the solute and solvent atom contributions, because
  # the data is read in a single column
  #for i in 1:S.input.nsave 
  #  for iat in 1:S.input.n
  #    traj.atoms[i].x[iat,1] = S.traj.atoms[i].x[iat]
  #    traj.atoms[i].x[iat,2] = S.traj.atoms[i].x[iat+S.input.n]
  #    traj.atoms[i].status[iat] = S.traj.atoms[i].status[iat]
  #  end
  #  traj.potential[i] = S.traj.potential[i]
  #  traj.kinetic[i] = S.traj.kinetic[i]
  #  traj.total[i] = S.traj.total[i]
  #  traj.time[i] = S.traj.time[i]
  #  traj.nenc[i] = S.traj.nenc[i]
  #  traj.U[i] = S.traj.U[i]
  #  traj.S[i] = S.traj.S[i]
  #  traj.D[i] = S.traj.D[i]
  #  traj.I[i] = S.traj.I[i]
  #end
  #input = S.input
  #S = nothing
  return R
end

