#
# Update the data with the data accumulated in a frame, for the non-self distributions
#
function update_counters_frame!(R :: Result, 
                                rdf_count_random_frame :: Vector{Float64}, 
                                md_count_random_frame :: Vector{Float64},
                                volume_frame :: Volume, solute :: Selection, solvent :: Selection,
                                n_solvent_in_bulk :: Float64)

  @. volume_frame.shell = 
        volume_frame.total * (rdf_count_random_frame/(R.options.n_random_samples*solvent.nmols))
  volume_frame.domain = sum(volume_frame.shell)
  # Volume of the bulk region
  if ! R.options.usecutoff
    volume_frame.bulk = volume_frame.total - volume_frame.domain
  else
    ibulk = setbin(R.dbulk+0.5*R.options.binstep,R.options.binstep)
    for i in ibulk:R.nbins
      volume_frame.bulk = volume_frame.bulk + volume_frame.shell[i]
    end
  end
  bulk_density_frame = (n_solvent_in_bulk/solute.nmols) / volume_frame.bulk 

  # The number of "random" molecules generated must be that to fill the complete
  # volume with the density of the bulk. Now that we know the density of the bulk
  # and the volume of the frame, we can adjust that random count
  fillup_factor = bulk_density_frame * volume_frame.total / solvent.nmols
  @. R.rdf_count_random = 
        R.rdf_count_random + rdf_count_random_frame*fillup_factor
  @. R.md_count_random = 
        R.md_count_random + md_count_random_frame*fillup_factor

  @. R.volume.shell = R.volume.shell + volume_frame.shell
  R.volume.bulk = R.volume.bulk + volume_frame.bulk
  R.volume.domain = R.volume.domain + volume_frame.domain
  R.density.solvent_bulk = R.density.solvent_bulk + bulk_density_frame

  return nothing
end

#
# Update counters with the data accumulated in a frame, for the self-distributions
#
#voltar: tem que acertar as normalizaçoes aqui
function update_counters_frame!(R :: Result, rdf_count_random_frame :: Vector{Float64},
                                volume_frame :: Volume, 
                                solvent :: Selection, 
                                npairs :: Int64, 
                                n_solvent_in_bulk :: Float64)

  @. R.rdf_count_random = R.rdf_count_random + rdf_count_random_frame
  @. volume_frame.shell = 
        volume_frame.total * (rdf_count_random_frame/(options.n_random_samples*solvent.nmols))
  volume_frame.domain = sum(volume_frame.shell)
  # Volume of the bulk region
  if ! R.options.usecutoff
    volume_frame.bulk = volume_frame.total - volume_frame.domain
  else
    ibulk = setbin(R.dbulk+0.5*R.options.binstep,R.options.binstep)
    for i in ibulk:R.nbins
      volume_frame.bulk = volume_frame.bulk + volume_frame.shell[i]
    end
  end

  @. R.volume.shell = R.volume.shell + volume_frame.shell
  R.volume.bulk = R.volume.bulk + volume_frame.bulk
  R.volume.domain = R.volume.domain + volume_frame.domain
  R.density.solvent_bulk = R.density.solvent_bulk + 
                           (solvent.nmols-1)*(n_solvent_in_bulk/npairs) / volume_frame.bulk

  return nothing
end

