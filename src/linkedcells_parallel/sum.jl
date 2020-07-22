#
# Sum the counts of two Results structures, adding the result to the first structure
# as in R1 = R1 + R2
#

function sum!( R1 :: Result, R2 :: Result )

  @. R1.md_count += R2.md_count
  @. R1.md_count_random += R2.md_count_random

  @. R1.solute_atom += R2.solute_atom
  @. R1.solvent_atom += R2.solvent_atom

  @. R1.rdf_count += R2.rdf_count
  @. R1.rdf_count_random += R2.rdf_count_random

  sum!(R1.density,R2.density)
  sum!(R1.volume,R2.volume)

end

function sum!( D1 :: Density, D2 :: Density )
  D1.solute += D2.solute
  D1.solvent += D2.solvent
  D1.solvent_bulk += D2.solvent_bulk
end

function sum!( V1 :: Volume, V2 :: Volume )
  V1.total += V2.total
  V1.bulk += V2.bulk
  V1.domain += V2.domain
  @. V1.shell += V2.shell
end
