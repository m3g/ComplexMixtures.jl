
using PDBTools 

function residue_contribution(R :: Result, pdb :: String; protein :: String = "protein") 

  protein_atoms = PDBTools.readPDB(atoms,protein)
  seq = PDBTools.getseq(protein_atoms)[1]

  nr = length(seq)
  nd = length(R.d)
  residue_contribution = zeros(nd,nr)

  for ir in 1:nr
    residue_atoms = selindex(protein,"resnum = $iresidue")
    for id in 1:nd
      residue_contribution[id,ir] = sum(R.solute_contrib[id,residue_atoms])
    end
  end

  return residue_contribution

end




