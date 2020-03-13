

solute = Solute( VMDselect("system.pdb","protein",vmd="/usr/bin/vmd"), 
                 nmols=1 )
solvent = Solvent( VMDselect("system.pdb","resname UREA",vmd="/usr/bin/vmd"), 
                   natomspermol=8 ) 





