"""
    merge(r::Vector{Result})

Merges the results of MDDF calculations obtained by running the same
analysis on multiple trajectories, or multiple parts of the same trajectory. It returns
a Result structure of the same type, with all the functions and counters representing averages
of the set provided weighted by the number of frames read in each Result set.

"""
function Base.merge(results::Vector{<:Result})
    cannot_merge = false
    nframes_read = 0
    for ir in eachindex(results)
        for file in results[ir].files
            nframes_read += file.nframes_read
        end
        for jr in ir+1:lastindex(results)
            if results[ir].nbins != results[jr].nbins
                println(
                    "ERROR: To merge Results, the number of bins of the histograms of both sets must be the same.",
                )
                cannot_merge = true
            end
            if !(results[ir].cutoff ≈ results[ir].cutoff)
                println(
                    "ERROR: To merge Results, cutoff distance of the of the histograms of both sets must be the same.",
                )
                cannot_merge = true
            end
            if (results[ir].solute != results[jr].solute) || (results[ir].solvent != results[jr].solvent)
                println(
                    "ERROR: To merge Results, the solute and solvent selections of both sets must be the same.",
                )
                cannot_merge = true
            end
        end
    end
    if cannot_merge
        throw(ArgumentError(" Incompatible set of results to merge. "))
    end

    # Total number of frames of all results, and total weight of frames
    nfiles = 0
    ntot_frames = 0
    tot_frame_weight = 0.0
    for result in results
        nfiles += length(result.files)
        ntot_frames += sum(file.nframes_read for file in result.files)
        tot_frame_weight += sum_frame_weights(result)
    end

    # Compute weight of each file, given the number frames. Also compute 
    # the weight of the data of each file, given the weights of the frames
    weights = zeros(nfiles)
    ifile = 0
    for result in results
        for file in result.files
            ifile += 1
            weights[ifile] = file.nframes_read / ntot_frames
        end
    end

    # Append all TrajectoryFileOptions to single vector
    files = TrajectoryFileOptions[]
    for result in results
        for file in result.files
            push!(files, file)
        end
    end

    # Initialize group counts
    solute_group_count = [zeros(results[1].nbins) for _ in 1:length(results[1].solute_group_count)]
    solvent_group_count = [zeros(results[1].nbins) for _ in 1:length(results[1].solvent_group_count)]
    solute_group_count_random = [zeros(results[1].nbins) for _ in 1:length(results[1].solute_group_count)]
    solvent_group_count_random = [zeros(results[1].nbins) for _ in 1:length(results[1].solvent_group_count)]

    # Structure for merged results
    R = Result(
        nbins=results[1].nbins,
        dbulk=results[1].dbulk,
        cutoff=results[1].cutoff,
        autocorrelation=results[1].autocorrelation,
        solute=results[1].solute,
        solvent=results[1].solvent,
        solute_group_count=solute_group_count,
        solvent_group_count=solvent_group_count,
        solute_group_count_random=solute_group_count_random,
        solvent_group_count_random=solvent_group_count_random,
        files=files,
        weights=weights,
    )

    # Average results weighting the data considering the weights of the frames of each data set
    warn = false
    @. R.d = results[1].d
    for ifile in eachindex(results)
        result = results[ifile]
        w = sum_frame_weights(result) / tot_frame_weight
        if !(w ≈ R.weights[ifile]) && !warn
            warn = true
            @warn begin
                """\n
                    Frame weights and file weights differ, because crustom frame weights were provided.

                """
            end _file = nothing _line = nothing
        end
        @. R.mddf += w * result.mddf
        @. R.kb += w * result.kb
        @. R.rdf += w * result.rdf
        @. R.kb_rdf += w * result.kb_rdf
        @. R.md_count += w * result.md_count
        @. R.md_count_random += w * result.md_count_random
        @. R.coordination_number += w * result.coordination_number
        @. R.coordination_number_random += w * result.coordination_number_random
        for i in eachindex(R.solute_group_count, result.solute_group_count)
            @. R.solute_group_count[i] += w * result.solute_group_count[i]
            @. R.solute_group_count_random[i] += w * result.solute_group_count_random[i]
        end
        for i in eachindex(R.solvent_group_count, result.solvent_group_count)
            @. R.solvent_group_count[i] += w * result.solvent_group_count[i]
            @. R.solvent_group_count_random[i] += w * result.solvent_group_count_random[i]
        end
        @. R.rdf_count += w * result.rdf_count
        @. R.rdf_count_random += w * result.rdf_count_random
        @. R.sum_rdf_count += w * result.sum_rdf_count
        @. R.sum_rdf_count_random += w * result.sum_rdf_count_random
        R.density.solute += w * result.density.solute
        R.density.solvent += w * result.density.solvent
        R.density.solvent_bulk += w * result.density.solvent_bulk
        R.volume.total += w * result.volume.total
        R.volume.bulk += w * result.volume.bulk
        R.volume.domain += w * result.volume.domain
        R.volume.shell += w * result.volume.shell
    end
    return R
end

@testitem "merge" begin
    using ComplexMixtures: mddf, merge
    using PDBTools: read_pdb, select, selindex
    using ComplexMixtures: data_dir

    # Test simple three-molecule system
    atoms = read_pdb("$data_dir/toy/cross.pdb")
    protein = AtomSelection(select(atoms, "protein and model 1"), nmols=1)
    water = AtomSelection(select(atoms, "resname WAT and model 1"), natomspermol=3)
    traj = Trajectory("$data_dir/toy/cross.pdb", protein, water, format="PDBTraj")

    options = Options(
        seed=321,
        StableRNG=true,
        nthreads=1,
        silent=true,
        n_random_samples=10^5,
        lastframe=1,
    )
    R1 = mddf(traj, options)

    options = Options(
        seed=321,
        StableRNG=true,
        nthreads=1,
        silent=true,
        n_random_samples=10^5,
        firstframe=2,
    )
    R2 = mddf(traj, options)

    R = merge([R1, R2])
    @test R.volume.total == 27000.0
    @test R.volume.domain ≈ R.volume.total - R.volume.bulk
    @test isapprox(R.volume.domain, (4π / 3) * R.dbulk^3; rtol=0.01)
    @test R.density.solute ≈ 1 / R.volume.total
    @test R.density.solvent ≈ 3 / R.volume.total
    @test R.density.solvent_bulk ≈ 2 / R.volume.bulk
    @test R.weights == [0.5, 0.5]

    # Test loading a saved merged file
    dir = mktempdir()
    save(R, "$dir/merged.json")
    R_save = load("$dir/merged.json")
    @test isapprox(R, R_save, debug=true)

    # Test merging files for which weights are provided for the frames
    R2 = mddf(traj, options, frame_weights=[0.0, 2.0])
    @test R.weights == [0.5, 0.5]
    @test length(R.files) == 2

    # Two-atom system
    at1 = AtomSelection([1], nmols=1)
    at2 = AtomSelection([2], nmols=1)
    traj = Trajectory("$data_dir/toy/self_monoatomic.pdb", at1, at2, format="PDBTraj")
    R1 = mddf(traj, Options(lastframe=1))
    @test sum(R1.md_count) == 1
    R2 = mddf(traj, Options(firstframe=2))
    @test sum(R2.md_count) == 0
    R = merge([R1, R2])
    @test R.weights == [0.5, 0.5]
    @test sum(R.md_count) == 0.5
    R1 = mddf(traj, Options(lastframe=1), frame_weights=[2.0])
    @test sum(R1.md_count) == 1
    R = merge([R1, R2])
    @test sum(R.md_count) == 2 / 3
    @test sum(sum.(R1.solute_group_count)) == 1
    @test sum(sum.(R1.solvent_group_count)) == 1
    @test sum(sum.(R2.solute_group_count)) == 0
    @test sum(sum.(R2.solvent_group_count)) == 0
    @test sum(sum.(R.solute_group_count)) == 2 / 3
    @test sum(sum.(R.solvent_group_count)) == 2 / 3

    # Test throwing merging incompatible results
    protein = AtomSelection(select(atoms, "protein and model 1"), nmols=1)
    traj = Trajectory("$data_dir/toy/cross.pdb", protein, water, format="PDBTraj")
    R1 = mddf(traj, options)
    protein = AtomSelection(
        select(atoms, "protein and model 1"), nmols=1,
        group_names=["acidic"],
        group_atom_indices=[selindex(atoms, "protein and acidic")]
    )
    traj = Trajectory("$data_dir/toy/cross.pdb", protein, water, format="PDBTraj")
    R2 = mddf(traj, options)


end
