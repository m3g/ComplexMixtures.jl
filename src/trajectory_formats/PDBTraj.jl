#
# Structure to contain PDB trajectories. Frames separated by "END", without any other field (remarks, etc)
#
# Must be mutable such that nframes, vector sizes, and additional parameters can be updated
# when the trajectory is first open
#

mutable struct PDBTraj

  # Mandatory data for things to work
  filename :: String
  iostream :: IOStream
  nframes :: Int64 

  # This Array (nframes,3) must be filled up with the size of the periodic cell
  sides :: Array{Float64}

  # Coordinates of the solute and solvent atoms in a frame (natoms,3) for each array:
  x_solute :: Array{Float64}  # (solute.n,3)
  x_solvent :: Array{Float64} # (solvent.n,3)

  # Additional data required for input/output functions
  natoms :: Int64

end

#
# Function open will set up the IO stream of the trajectory, fillup the 
# number of frames field and additional parameters if required (in this case,
# if the trajectory file has or not pbc information). 
#
# The IO stream must be returned in a position such that the "nextframe" function
# will be able to read the first frame of the trajectory
#

function open( trajectory :: PDBTraj )

  pdbfile = open(trajectory.filename,"r")
  
  # Get the number of frames and number of atoms, check for adecuacy of the format
  nframes = 0
  natoms = 0
  for line in eachline(pdbfile)
    if line(1:3) == "END"
      nframes = nframes + 1
    elseif nframes = 0 && (line(1:4) == "ATOM" || line(1:6) == "HETATM")
      natoms = natoms + 1
    else
      error(" The PDB file containing the trajectory must have only ATOM, HETATM and END fields. ")
    end
  end
  close(pdbfile)

  trajectory.nframes = nframes
  trajectory.natoms = natoms

  # Open again the iostream, so that nextrame can read the first frame
  trajectory.iostream = open(trajectory.filename,"r")
  
end

#
# Function that reads the coordinates of the solute and solvent atoms from
# the next frame of the trajectory file 
#
# The function modifies sides, x_solute and x_solvent within the trajectory structure.
# Having these vectors inside the trajectory structure avoids having to allocate
# them everytime a new frame is read
#

function nextframe!( trajectory:: PDBTraj, solute :: Solute, solvent :: Solvent )

  iatom = 0
  line = "START"
  x = Array{Float64}(undef,natoms,3)
  while line(1:3) != "END"
    line = readline(trajectory.iostream)
    iatom = iatom + 1
    x[iatom,1] = parse(Float64,record[31:38])
    x[iatom,2] = parse(Float64,record[39:46])
    x[iatom,3] = parse(Float64,record[47:54])
  end

  # Save coordinates of solute and solvent in trajectory arrays
  for i in 1:solute.n
    trajectory.x_solute[i,1] = x[solute.index[i],1]
    trajectory.x_solute[i,2] = x[solute.index[i],2]
    trajectory.x_solute[i,3] = x[solute.index[i],3]
  end
  for i in 1:solvent.n
    trajectory.x_solvent[i,1] = x[solvent.index[i],1]
    trajectory.x_solvent[i,2] = x[solvent.index[i],2]
    trajectory.x_solvent[i,3] = x[solvent.index[i],3]
  end

end

#
# Function that closes the IO Stream of the trajectory
#

function close( trajectory :: PDBTraj )
  close(trajectory.iostream)
end

