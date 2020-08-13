#
# Function to print a resume of the results on the screen after running
#

function print_results( results :: Results )

  println()
  println(@printf("%s %12.5f","  Solvent density in simulation box (sites/A^3): ", simdensity))
  println(@printf("%s %12.5f","  Estimated bulk solvent density (sites/A^3): ", bulkdensity))
  println()                   
  println(@printf("%s %12.5f","  Molar volume of solvent in simulation box (cc/mol): ", convert/simdensity))
  println(@printf("%s %12.5f","  Molar volume of solvent in bulk (cc/mol): ", convert/bulkdensity))
  println()                   
  println(@printf("%s %12.5f","  Density scaling factor for numerical integration: ", density_fix))

  println()
  println(@printf("%s %12.5f","  Solute partial volume (cc/mol): ", solutevolume))

  return nothing
end
 

