#=

Structure to carry temporary arrays needed

We need temporary buffers for reading the coordinates, a buffer to store
a temporary array of solvent coordinates (whose length is dependent on 
autocorrelation or not), and the indices of the solute molecules that will
be used as reference states for ideal gas distributions

=#
@kwdef struct Buffer
    solute_read::Vector{SVector{3,Float64}}
    solvent_read::Vector{SVector{3,Float64}}
    ref_solutes::Vector{Int}
    list::Vector{MinimumDistance}
    indices_in_bulk::Vector{Int}
end
function Buffer(trajectory::Trajectory, R::Result)
    return Buffer(
        solute_read=similar(trajectory.x_solute),
        solvent_read=similar(trajectory.x_solvent),
        ref_solutes=zeros(Int, R.files[1].options.n_random_samples),
        list=fill(zero(MinimumDistance), R.solvent.nmols),
        indices_in_bulk=fill(0, R.solvent.nmols),
    )
end

@testitem "Buffer" begin
    using ComplexMixtures
    using PDBTools
    using ComplexMixtures.Testing
    atoms = readPDB(Testing.pdbfile)
    options = Options(stride=5, seed=321, StableRNG=true, nthreads=1, silent=true)
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    traj = Trajectory("$(Testing.data_dir)/NAMD/trajectory.dcd", protein, tmao)
    R = Result(traj, options)
    b0 = ComplexMixtures.Buffer(
        solute_read=similar(traj.x_solute),
        solvent_read=similar(traj.x_solvent),
        ref_solutes=zeros(Int, R.files[1].options.n_random_samples),
        list=fill(zero(ComplexMixtures.MinimumDistance), R.solvent.nmols),
        indices_in_bulk=fill(0, R.solvent.nmols),
    )
    b1 = ComplexMixtures.Buffer(traj, R)
    for field in fieldnames(ComplexMixtures.Buffer)
        @test typeof(getfield(b0, field)) == typeof(getfield(b1, field))
        @test length(getfield(b0, field)) == length(getfield(b1, field))
    end
end

#
# Defines is a molecule is a bulk molecule
#
inbulk(md::MinimumDistance, options::Options) =
    options.usecutoff ? (md.within_cutoff && md.d > options.dbulk) :
    !md.within_cutoff

#=
    randomize_solvent!(system, buff, n_solvent_in_bulk, options, RNG)

Generate a random solvent distribution from the bulk molecules of a solvent

=#
function randomize_solvent!(
    system::AbstractParticleSystem,
    buff::Buffer,
    n_solvent_in_bulk::Integer,
    R::Result,
    RNG,
)
    for isolvent = 1:R.solvent.nmols
        # Choose randomly one molecule from the bulk, if there are bulk molecules
        if n_solvent_in_bulk > 0
            jmol = buff.indices_in_bulk[rand(RNG, 1:n_solvent_in_bulk)]
        else
            jmol = rand(RNG, 1:R.solvent.nmols)
        end
        # Pick coordinates of the molecule to be randomly moved
        y_new = viewmol(isolvent, system.ypositions, R.solvent)
        # Copy the coordinates of the random solvent molecule chosen
        y_new .= viewmol(jmol, buff.solvent_read, R.solvent)
        # Randomize rotations and translation for this molecule 
        random_move!(y_new, R.files[1].irefatom, system, RNG)
    end
end

# 
# Read trajectory until next frame to be read. In a function to try to 
# improve the GC behavior, because there is a leakage in memory in the
# reading of the coordinates by Chemfiles
#
function goto_nextframe!(iframe, R, trajectory, to_compute_frames, options)
    frame_weight = 1.0
    compute = false
    while iframe < R.files[1].lastframe_read && !compute
        nextframe!(trajectory)
        iframe += 1
        frame_weight = isempty(R.files[1].frame_weights) ? 1.0 : R.files[1].frame_weights[iframe]
        if iframe in to_compute_frames && !iszero(frame_weight)
            compute = true
            # Run GC if memory is getting full: this is related to issues with Chemfiles reading scheme
            if options.GC && (Sys.free_memory() / Sys.total_memory() < options.GC_threshold)
                GC.gc()
            end
        end
    end
    return iframe, frame_weight, compute
end

"""     
    mddf(
        trajectory_file::String, 
        solute::AtomSelection,
        solvent::AtomSelection, # optional: if omitted, an auto-correlation will be computed
        options::Options=Options(); # optional: if omitted, default options will be used
        # optional keywords:
        trajectory_format="", 
        frame_weights = Float64[], 
        coordination_number_only = false, 
        low_memory = false, 
        chemfiles = false,
    )

Computes the minimum-distance distribution function, atomic contributions, and 
KB integrals, given the trajectory file name, and the definition of the solute and
solvent groups of atoms as `AtomSelection` data structures. If the `solvent` parameter
is omitted, a self-correlation will be computed.

This is the main function of the `ComplexMixtures` package. 

The `options` parameter is optional. If not set, the default `Options()` structure will be used.

### Optional execution keywords

- `trajectory_format` is a string that defines the format of the trajectory file. If not provided, the format
   will be guessed from the file extension, and the `Chemfiles` library will be used to read the trajectory. 

- `frame_weights` is an array of weights for each frame of the trajectory. If this is provided, the MDDF will be computed as a weighted
  average of the MDDFs of each frame. This can be used to study the solvation dependency in perturbed 
  ensembles.

- `coordination_number_only` is a boolean that, if set to `true`, will compute only the
  site-counts and coordination numbers of the solvent molecules around the solute, and not the MDDF, RDF, or KB integrals. 
  This is useful when the normalization of the distribution is not possible or needed, for instance when
  the bulk solution is not properly defined. The computation is much faster in this case. The call to `mddf` with this
  option is equivalent to direct the call to `coordination_number`.

- `low_memory` can be set to `true` to reduce the memory requirements of the computation. This will
  parallelize the computation of the minimum distances at a lower level, reducing the memory
  requirements at the expense of some performance.

- `chemfiles` is a boolean that, if set to `true`, will use try to use the `Chemfiles` library to read the trajectory file. 
  independently of the file extension. By default, PDB and DCD trajectories are read by specific readers.

!!! compat
    The current call signature was introduced in version 2.9.0. The previous call signature with the 
    previous construction of the `Trajectory` object, is still available, but it is no longer recommended 
    because it more prone to user-level solute and solvent configuration errors. 

### Examples

```julia-repl
julia> using ComplexMixtures, PDBTools

julia> using ComplexMixtures.Testing: data_dir

julia> atoms = readPDB(joinpath(data_dir,"NAMD/structure.pdb"));

julia> solute = AtomSelection(select(atoms, "protein"), nmols=1);

julia> solvent = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14);

julia> options = Options(lastframe=10, bulk_range=(10.0, 15.0));

julia> trajectory_file = joinpath(data_dir,"NAMD/trajectory.dcd");

julia> results = mddf(trajectory_file, solute, solvent, options);
```

"""
function mddf end

function mddf(
    trajectory_file::String,
    solute::AtomSelection,
    solvent::AtomSelection,
    options::Options=Options();
    trajectory_format::String="",
    chemfiles::Bool=false,
    kargs...)
    trajectory = Trajectory(trajectory_file, solute, solvent; format=trajectory_format, chemfiles)
    return mddf(trajectory, options; kargs...)
end

function mddf(
    trajectory_file::String,
    solute_and_solvent::AtomSelection,
    options::Options=Options();
    trajectory_format::String="",
    chemfiles::Bool=false,
    kargs...
)
    trajectory = Trajectory(trajectory_file, solute_and_solvent; format=trajectory_format, chemfiles)
    return mddf(trajectory, options; kargs...)
end

function mddf(
    trajectory::Trajectory, options::Options=Options();
    frame_weights=Float64[],
    coordination_number_only=false,
    low_memory=false,
)

    options.silent || println(bars)
    options.silent || println("Initializing data structures ...")

    # Set random number generator
    RNG = init_random(options)

    # Compute some meta-data from the trajectory file and the options,
    # to allow parallel creation of the Result, and ParticleSystem
    # data structures
    options.silent || println("Evaluate trajectory length and go to first frame ...")
    trajectory_data = TrajectoryMetaData(trajectory, options)

    # Structure in which the final results will be stored
    R = Result(trajectory, options; trajectory_data, frame_weights)

    # Open the trajectory stream and go to first frame
    opentraj!(trajectory)
    firstframe!(trajectory)

    # Skip initial frames if desired
    progress = Progress(options.firstframe; dt=1, enabled=!options.silent)
    for _ in 1:options.firstframe-1
        nextframe!(trajectory)
        if options.GC && (Sys.free_memory() / Sys.total_memory() < options.GC_threshold)
            GC.gc()
        end
        next!(progress)
    end
    iframe = options.firstframe - 1

    # Define how the parallelization will be performed, according to the memory
    # requirements of the computation
    nchunks, parallel_cl, nbatches_cl, nthreads = parallel_setup(options, R, low_memory)

    # Print some information about this run
    if !options.silent
        title(R, trajectory.solute, trajectory.solvent, nthreads)
        if low_memory
            println("Running with low-memory option:")
            println("  - Parallel CellListMap can be used.")
            println("  - Number of parallel minimum-distance computations: $nchunks")
            println("  - Each minimum-distance computation will use $(nbatches_cl[2]) threads.")
        end
    end

    # Frames to be read and frames for which the MDDF will be computed
    to_read_frames = options.firstframe:R.files[1].lastframe_read
    to_compute_frames = options.firstframe:options.stride:R.files[1].lastframe_read
    progress = Progress(
        count(i -> isempty(frame_weights) ? true : !iszero(frame_weights[i]), to_compute_frames); 
        dt=1, enabled=!options.silent
    )

    # Loop over the trajectory
    read_lock = ReentrantLock()
    sum_results_lock = ReentrantLock()
    @sync for frame_range in ChunkSplitters.chunks(to_read_frames; n=nchunks)
        Threads.@spawn begin
            # Local data structures for this chunk
            system_chunk = ParticleSystem(
                trajectory, trajectory_data.unitcell, options, parallel_cl, nbatches_cl
            )
            buff_chunk = Buffer(trajectory, R)
            r_chunk = Result(trajectory, options; trajectory_data, frame_weights)
            # Reset the number of frames read by each chunk
            for _ in frame_range
                # interrupt handler
                if isfile("stop_complexmixtures") 
                    @info "From thread id: $(Threads.threadid()): stop_complexmixtures file found. Exiting."
                    break
                end
                local compute, frame_weight
                # Read frame coordinates
                @lock read_lock begin
                    iframe, frame_weight, compute = goto_nextframe!(iframe, R, trajectory, to_compute_frames, options)
                    if compute
                        # Read frame for computing 
                        # The solute coordinates must be read in intermediate arrays, because the 
                        # solute molecules will be considered one at a time in the computation of the
                        # minimum distances
                        @. buff_chunk.solute_read = trajectory.x_solute
                        @. buff_chunk.solvent_read = trajectory.x_solvent
                        # Read weight of this frame
                        unitcell = convert_unitcell(trajectory_data.unitcell, getunitcell(trajectory))
                        update_unitcell!(system_chunk, unitcell)
                        # Display progress bar
                        next!(progress)
                    end
                end # release reading lock
                #
                # Perform MDDF computation
                #
                if compute
                    # Compute distances in this frame and update results
                    if !coordination_number_only
                        mddf_frame!(r_chunk, system_chunk, buff_chunk, options, frame_weight, RNG)
                    else
                        coordination_number_frame!(r_chunk, system_chunk, buff_chunk, frame_weight)
                    end
                end # compute if
            end # frame range for this chunk
            # Sum the results of this chunk into the total results data structure
            @lock sum_results_lock sum!(R, r_chunk)
        end # spawn
    end # chunk loop
    closetraj!(trajectory)

    # Setup the final data structure with final values averaged over the number of frames,
    # sampling, etc, and computes final distributions and integrals
    finalresults!(R, options; coordination_number_only)
    options.silent || println(bars)

    return R
end

# Compute cell volume from unitcell matrix
cell_volume(system::AbstractParticleSystem) =
    @views dot(cross(system.unitcell[:, 1], system.unitcell[:, 2]), system.unitcell[:, 3])
update_volume!(r_chunk, system, frame_weight) = r_chunk.volume.total += frame_weight * cell_volume(system)

#=
    mddf_frame!(r_chunk::Result, system::AbstractParticleSystem, buff::Buffer, options::Options, frame_weight, RNG)

Computes the MDDF for a single frame. Modifies the data in the `r_chunk` (type `Result`) structure, through
the `update_volume!` and `update_counters!` functions.

=#
function mddf_frame!(
    r_chunk::Result,
    system::AbstractParticleSystem,
    buff::Buffer,
    options::Options,
    frame_weight,
    RNG,
)

    # Sum up the volume of this frame
    update_volume!(r_chunk, system, frame_weight)

    # Random set of solute molecules to use as reference for the ideal gas distributions
    for i in eachindex(buff.ref_solutes)
        buff.ref_solutes[i] = rand(RNG, 1:r_chunk.solute.nmols)
    end

    #
    # Compute the MDDFs for each solute molecule
    #
    update_lists = true
    system.ypositions .= buff.solvent_read
    for isolute = 1:r_chunk.solute.nmols

        # We need to do this one solute molecule at a time to avoid exploding the memory requirements
        system.xpositions .= viewmol(isolute, buff.solute_read, r_chunk.solute)

        # Compute minimum distances of the molecules to the solute (updates system.list, and returns it)
        # The cell lists will be recomputed for the first solute, or if a random distribution was computed
        # for the previous solute
        minimum_distances!(system, r_chunk, isolute; update_lists=update_lists)
        update_lists = false # avoid recomputation of the cell lists in the next iteration

        # For each solute molecule, update the counters (this is highly suboptimal, because
        # within update_counters there are loops over solvent molecules, in such a way that
        # this will loop with cost nsolute*nsolvent. However, I cannot see an easy solution 
        # at this point with acceptable memory requirements
        update_counters!(r_chunk, system, frame_weight)

        # If this molecule was chosen as a reference molecule for the ideal gas distribution, compute it
        # (as many times as needed, as the reference molecules may be repeated - particularly because
        # there may be only one solute molecule, in which case all distributions will be created for
        # the same solute molecule).
        nrand = count(==(isolute), buff.ref_solutes)
        if nrand > 0
            update_lists = true
            # Save coordinates and list
            buff.list .= system.list
            # Annotate the indices of the molecules that are in the bulk solution
            n_solvent_in_bulk = 0
            for i in eachindex(system.list)
                r_chunk.autocorrelation && i == isolute && continue
                if inbulk(system.list[i], options)
                    n_solvent_in_bulk += 1
                    buff.indices_in_bulk[n_solvent_in_bulk] = i
                end
            end
            for _ = 1:nrand
                randomize_solvent!(system, buff, n_solvent_in_bulk, r_chunk, RNG)
                minimum_distances!(system, r_chunk, isolute; update_lists=update_lists)
                update_counters!(r_chunk, system, frame_weight, Val(:random))
            end
            # Restore positions and minimum distance list of system structure
            system.list .= buff.list
            system.ypositions .= buff.solvent_read
        end

    end # loop over solute molecules

    return r_chunk
end

#=
    coordination_number_frame!(r_chunk::Result, system::AbstractParticleSystem, buff::Buffer, frame_weight)

Computes the coordination number for a single frame. Modifies the data in the `r_chunk` (type `Result`) structure, through
the `update_volume!` and `update_counters!` functions.

=#
function coordination_number_frame!(
    r_chunk::Result,
    system::AbstractParticleSystem,
    buff::Buffer,
    frame_weight::Float64,
)

    # Sum up the volume of this frame
    update_volume!(r_chunk, system, frame_weight)

    #
    # Compute the MDDFs for each solute molecule
    #
    update_lists = true
    system.ypositions .= buff.solvent_read
    for isolute = 1:r_chunk.solute.nmols

        # We need to do this one solute molecule at a time to avoid exploding the memory requirements
        system.xpositions .= viewmol(isolute, buff.solute_read, r_chunk.solute)

        # Compute minimum distances of the molecules to the solute (updates system.list, and returns it)
        # The cell lists will be recomputed for the first solute, or if a random distribution was computed
        # for the previous solute
        minimum_distances!(system, r_chunk, isolute; update_lists=update_lists)
        update_lists = false # avoid recomputation of the cell lists in the next iteration

        # For each solute molecule, update the counters (this is highly suboptimal, because
        # within update_counters there are loops over solvent molecules, in such a way that
        # this will loop with cost nsolute*nsolvent. However, I cannot see an easy solution 
        # at this point with acceptable memory requirements
        update_counters!(r_chunk, system, frame_weight)

    end # loop over solute molecules

    return r_chunk
end


"""     
    coordination_number(
        trajectory_file::String, 
        solute::AtomSelection,
        solvent::AtomSelection, # optional: if omitted, an auto-correlation will be computed
        options::Options=Options(); # optional: if omitted, default options will be used
        # optional keywords:
        trajectory_format="", 
        frame_weights = Float64[], 
        low_memory = false, 
    )

Computes the minimum-distance coordination number and atomic contributions,
given the trajectory file name, and the definition of the solute and
solvent groups of atoms as `AtomSelection` data structures. If the `solvent` parameter
is omitted, a self-coordination number will be computed.

The output is a `Result` structure, which contains the data as the result of a call to `mddf`, except
that all counters which require normalization of the distribution will be zero. In summary, this result
data structure can be used to compute the coordination numbers, but not the MDDF, RDF, or KB integrals.

The `options` parameter is optional. If not set, the default `Options()` structure will be used.

### Optional execution keywords

- `trajectory_format` is a string that defines the format of the trajectory file. If not provided, the format
   will be guessed from the file extension, and the `Chemfiles` library will be used to read the trajectory. 

- `frame_weights` is an array of weights for each frame of the trajectory. If this is provided, the MDDF will be computed as a weighted
  average of the MDDFs of each frame. This can be used to study the solvation dependency in perturbed 
  ensembles.

- `low_memory` can be set to `true` to reduce the memory requirements of the computation. This will
  parallelize the computation of the minimum distances at a lower level, reducing the memory
  requirements at the expense of some performance.

- `chemfiles` is a boolean that, if set to `true`, will use try to use the `Chemfiles` library to read the trajectory file. 
  independently of the file extension. By default, PDB and DCD trajectories are read by specific readers.

!!! compat
    The current call signature was introduced in version 2.9.0. The previous call signature with the 
    previous construction of the `Trajectory` object, is still available, but it is no longer recommended 
    because it more prone to user-level solute and solvent configuration errors. 

### Examples

```julia-repl
julia> using ComplexMixtures, PDBTools

julia> using ComplexMixtures.Testing: data_dir

julia> atoms = readPDB(joinpath(data_dir,"NAMD/structure.pdb"));

julia> solute = AtomSelection(select(atoms, "protein"), nmols=1);

julia> solvent = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14);

julia> options = Options(lastframe=10, bulk_range=(10.0, 15.0));

julia> trajectory_file = joinpath(data_dir,"NAMD/trajectory.dcd");

julia> results = coordination_number(trajectory_file, solute, solvent, options);
```

"""
function coordination_number(
    trajectory_file::String,
    solute::AtomSelection,
    solvent::AtomSelection,
    options::Options=Options();
    trajectory_format::String="",
    chemfiles::Bool=false,
    kargs...)
    _coordination_number_call_error(; kargs...)
    trajectory = Trajectory(trajectory_file, solute, solvent; format=trajectory_format, chemfiles)
    return mddf(trajectory, options; coordination_number_only=true, kargs...)
end

function coordination_number(
    trajectory_file::String,
    solute_and_solvent::AtomSelection,
    options::Options=Options();
    trajectory_format::String="",
    chemfiles::Bool=false,
    kargs...
)
    _coordination_number_call_error(; kargs...)
    trajectory = Trajectory(trajectory_file, solute_and_solvent; format=trajectory_format, chemfiles)
    return mddf(trajectory, options; coordination_number_only=true, kargs...)
end

function coordination_number(traj::Trajectory, options::Options=Options(); kargs...)
    _coordination_number_call_error(; kargs...)
    return mddf(traj, options; coordination_number_only=true, kargs...)
end

function _coordination_number_call_error(; kargs...)
    if haskey(kargs, :coordination_number_only)
        throw(ArgumentError("""\n
            The keyword argument `coordination_number_only` is not allowed in the `coordination_number` function.
            It is, by definition, set to `true` in this function. 
            
        """))
    end
    return nothing
end

@testitem "mddf - toy system" begin
    using ComplexMixtures
    using PDBTools: readPDB, select
    using ComplexMixtures.Testing: data_dir

    # Test simple three-molecule system: cross correlation
    atoms = readPDB("$data_dir/toy/cross.pdb")
    protein = AtomSelection(select(atoms, "protein and model 1"), nmols=1)
    water = AtomSelection(select(atoms, "resname WAT and model 1"), natomspermol=3)
    trajectory_file = "$data_dir/toy/cross.pdb"
    trajectory_format = "PDBTraj"
    for nthreads in [1, 2], lastframe in [1, 2], low_memory in [true, false]
        options = Options(;
            seed=321,
            StableRNG=true,
            nthreads,
            silent=true,
            n_random_samples=10^5,
            lastframe,
        )
        R = mddf(trajectory_file, protein, water, options; trajectory_format, low_memory)
        @test R.volume.total == 27000.0
        @test R.volume.domain ≈ R.volume.total - R.volume.bulk
        @test isapprox(R.volume.domain, (4π / 3) * R.dbulk^3; rtol=0.01)
        @test R.density.solute ≈ 1 / R.volume.total
        @test R.density.solvent ≈ 3 / R.volume.total
        @test R.density.solvent_bulk ≈ 2 / R.volume.bulk
        @test sum(R.md_count) ≈ 1
        @test sum(R.coordination_number) ≈ 51
        C = coordination_number(trajectory_file, protein, water, options; trajectory_format, low_memory)
        @test C.volume.total == R.volume.total
        @test C.volume.domain ≈ 0.0
        @test C.density.solute ≈ 1 / C.volume.total
        @test C.density.solvent ≈ 3 / C.volume.total
        @test C.density.solvent_bulk == 0.0
        @test C.md_count == R.md_count
        @test coordination_number(C) == coordination_number(R)
    end

    # Test wrong frame_weights input
    @test_throws ArgumentError mddf(trajectory_file, protein, water, Options(); frame_weights=[1.0], trajectory_format)

    # Self correlation
    atoms = readPDB("$data_dir/toy/self_monoatomic.pdb")
    atom = AtomSelection(select(atoms, "resname WAT and model 1"), natomspermol=1)
    trajectory_file = "$data_dir/toy/self_monoatomic.pdb"
    trajectory_format = "PDBTraj"
    # without atoms in the bulk
    for nthreads in [1, 2], low_memory in [false, true]
        options = Options(;
            seed=321,
            StableRNG=true,
            nthreads,
            silent=true,
            n_random_samples=10^5,
            lastframe=1,
        )
        R = mddf(trajectory_file, atom, options; low_memory, trajectory_format)
        @test R.volume.total == 27000.0
        @test R.volume.domain ≈ R.volume.total - R.volume.bulk
        @test isapprox(R.volume.domain, (4π / 3) * R.dbulk^3; rtol=0.1)
        @test R.density.solute ≈ 2 / R.volume.total
        @test R.density.solvent ≈ 2 / R.volume.total
        @test sum(R.md_count) ≈ 1.0
        R = coordination_number(trajectory_file, atom, options; low_memory, trajectory_format)
        @test sum(R.md_count) ≈ 1.0

    end

    # Test interruption
    current_dir = pwd()
    test_dir = tempname()
    mkdir(test_dir)
    cd(test_dir)
    touch("stop_complexmixtures")
    R_interrupt = mddf(trajectory_file, atom)
    @test sum(R_interrupt.md_count) == 0.0
    cd(current_dir)

    #
    # Test varying frame weights
    #
    # Read only first frame
    for low_memory in [false, true]
        local trajectory_file
        options = Options(seed=321, StableRNG=true, nthreads=1, silent=true, lastframe=1)
        trajectory_file = "$data_dir/toy/self_monoatomic.pdb"
        R1 = mddf(trajectory_file, atom, options; low_memory, trajectory_format)
        R2 = mddf(trajectory_file, atom, options; frame_weights=[1.0, 0.0], low_memory, trajectory_format)
        @test R1.md_count == R2.md_count
        @test R1.rdf_count == R2.rdf_count
        # Read only last frame
        options = Options(seed=321, StableRNG=true, nthreads=1, silent=true, firstframe=2)
        R1 = mddf(trajectory_file, atom, options; low_memory, trajectory_format)
        options = Options(seed=321, StableRNG=true, nthreads=1, silent=true)
        R2 = mddf(trajectory_file, atom, options; frame_weights=[0.0, 1.0], low_memory, trajectory_format)
        @test R1.md_count == R2.md_count
        @test R1.rdf_count == R2.rdf_count
        # Use equal weights
        R1 = mddf(trajectory_file, atom, options; low_memory, trajectory_format)
        R2 = mddf(trajectory_file, atom, options; frame_weights=[0.3, 0.3], low_memory, trajectory_format)
        @test R1.md_count == R2.md_count
        @test R1.rdf_count == R2.rdf_count
        # Check with the duplicated-first-frame trajectory 
        trajectory_file = "$data_dir/toy/self_monoatomic_duplicated_first_frame.pdb"
        options = Options(seed=321, StableRNG=true, nthreads=1, silent=true, firstframe=1, lastframe=3)
        R1 = mddf(trajectory_file, atom, options; low_memory, trajectory_format)
        trajectory_file = "$data_dir/toy/self_monoatomic.pdb"
        options = Options(seed=321, StableRNG=true, nthreads=1, silent=true)
        R2 = mddf(trajectory_file, atom, options; frame_weights=[2.0, 1.0], low_memory, trajectory_format)
        @test R1.md_count == R2.md_count
        @test R1.rdf_count == R2.rdf_count
        # Test some different weights
        R2 = mddf(trajectory_file, atom, Options(silent=true); frame_weights=[2.0, 1.0], low_memory, trajectory_format)
        @test sum(R2.md_count) ≈ 2 / 3
        R2 = mddf(trajectory_file, atom, Options(silent=true); frame_weights=[1.0, 2.0], low_memory, trajectory_format)
        @test sum(R2.md_count) ≈ 1 / 3

        # only with atoms in the bulk
        options = Options(
            seed=321,
            StableRNG=true,
            nthreads=1,
            silent=true,
            n_random_samples=10^5,
            firstframe=2,
        )
        R = mddf(trajectory_file, atom, options; low_memory, trajectory_format)
        @test R.volume.total == 27000.0
        @test R.volume.domain ≈ R.volume.total - R.volume.bulk
        @test isapprox(R.volume.domain, (4π / 3) * R.dbulk^3; rtol=0.1)
        @test R.density.solute ≈ 2 / R.volume.total
        @test R.density.solvent ≈ 2 / R.volume.total
        @test R.density.solvent_bulk ≈ 1 / R.volume.bulk

        # with both frames
        options = Options(
            seed=321,
            StableRNG=true,
            nthreads=1,
            silent=true,
            n_random_samples=10^5,
        )
        R = mddf(trajectory_file, atom, options; low_memory, trajectory_format)
        @test R.volume.total == 27000.0
        @test R.volume.domain ≈ R.volume.total - R.volume.bulk
        @test isapprox(R.volume.domain, (4π / 3) * R.dbulk^3; rtol=0.1)
        @test R.density.solute ≈ 2 / R.volume.total
        @test R.density.solvent ≈ 2 / R.volume.total
        @test R.density.solvent_bulk ≈ 0.5 / R.volume.bulk
    end
end

@testitem "mddf - available methods" begin
    using ComplexMixtures
    using PDBTools: readPDB, select
    using ComplexMixtures.Testing: data_dir
    atoms = readPDB("$data_dir/toy/cross.pdb")
    protein = AtomSelection(select(atoms, "protein and model 1"), nmols=1)
    water = AtomSelection(select(atoms, "resname WAT and model 1"), natomspermol=3)
    trajectory_file = "$data_dir/toy/cross.pdb"
    
    atom = AtomSelection(select(atoms, "resname WAT and model 1"), natomspermol=1)
    t = Trajectory(trajectory_file, atom; format="PDBTraj")

    options = Options()
    trajectory_file = "$data_dir/toy/self_monoatomic.pdb"
    atom = AtomSelection(select(atoms, "resname WAT and model 1"), natomspermol=1)
    r1 = mddf(trajectory_file, atom, atom)
    r2 = mddf(trajectory_file, atom, atom, Options())
    @test coordination_number(r1) ≈ coordination_number(r2)
    r1 = mddf(trajectory_file, atom)
    r2 = mddf(trajectory_file, atom, Options())
    @test coordination_number(r1) ≈ coordination_number(r2)
    r1 = coordination_number(trajectory_file, atom, atom)
    r2 = coordination_number(trajectory_file, atom, atom, Options())
    @test coordination_number(r1) ≈ coordination_number(r2)
    r1 = coordination_number(trajectory_file, atom)
    r2 = coordination_number(trajectory_file, atom, Options())
    @test coordination_number(r1) ≈ coordination_number(r2)
end

@testitem "mddf - real system" begin
    using ComplexMixtures: mddf, Trajectory, Options, AtomSelection
    using PDBTools: readPDB, select
    using ComplexMixtures.Testing: data_dir, pdbfile

    atoms = readPDB(pdbfile)
    protein = AtomSelection(select(atoms, "protein"), nmols=1)
    tmao = AtomSelection(select(atoms, "resname TMAO"), natomspermol=14)
    trajectory_file = "$data_dir/NAMD/trajectory.dcd"

    # Throw error in incorrect call to coordination_number
    @test_throws ArgumentError coordination_number(trajectory_file, protein, tmao, Options(); coordination_number_only=true)
    @test_throws ArgumentError coordination_number(trajectory_file, tmao, Options(); coordination_number_only=true)

    # Throw insufficent memory error
    @test_throws ErrorException mddf(trajectory_file, protein, tmao, Options(nthreads=10^10))

    for nthreads in [1, 2], low_memory in [true, false]
        options = Options(;
            seed=1,
            stride=1,
            StableRNG=true,
            nthreads,
            n_random_samples=100,
            silent=true
        )
        # Test actual system: cross correlation
        R = mddf("$data_dir/NAMD/trajectory.dcd", protein, tmao, options; low_memory)
        # Deterministic
        @test R.volume.total ≈ 603078.4438609097
        @test sum(R.md_count) ≈ 25.250000000000007
        @test sum(R.rdf_count) ≈ 19.550000000000004
        @test R.density.solute ≈ 1.6581590839128614e-6
        @test R.density.solvent ≈ 0.00030012679418822794
        # Dependent on the random number seed
        @test R.volume.domain ≈ 74560.15401932324 rtol = 0.1
        @test R.volume.bulk ≈ 528518.2898415865 rtol = 0.1
        @test R.density.solvent_bulk ≈ 0.0003054766563488117 rtol = 0.1
        @test sum(R.mddf) ≈ 527.670164155438 rtol = 0.1
        @test sum(R.rdf) ≈ 444.71561185073836 rtol = 0.1
        @test R.kb[end] ≈ -5775.215792514756 rtol = 0.5
        @test R.kb_rdf[end] ≈ -6360.471034166915 rtol = 0.5

        # Self correlation
        R = mddf("$data_dir/NAMD/trajectory.dcd", tmao, options; low_memory)
        # Deterministic
        @test R.volume.total ≈ 603078.4438609097
        @test sum(R.md_count) ≈ 2.8939226519337016
        @test sum(R.rdf_count) ≈ 1.858839779005525
        @test R.density.solute ≈ 0.00030012679418822794
        @test R.density.solvent ≈ 0.00030012679418822794
        # Dependent on the random number seed
        if nthreads == 1
            @test R.volume.domain ≈ 6958.855154995052 rtol = 0.1
            @test R.volume.bulk ≈ 596277.0591884783 rtol = 0.1
            @test R.density.solvent_bulk ≈ 0.00029875568324470034 rtol = 0.1
            @test sum(R.mddf) ≈ 434.7773107740875 rtol = 0.1
            @test sum(R.rdf) ≈ 345.9169130954568 rtol = 0.1
            @test R.kb[end] ≈ -476 rtol = 0.5
            @test R.kb_rdf[end] ≈ -421 rtol = 0.5
        end

        # Test varying frame weights: the trajectory below has 3 frames
        # extracted from NAMD/trajectory.dcd, and the 2 first frames are the
        # first frame duplicated.
        R1 = mddf("$data_dir/NAMD/traj_duplicated_first_frame.dcd", tmao, Options(); low_memory)
        R2 = mddf("$data_dir/NAMD/trajectory.dcd", tmao, Options(lastframe=2); frame_weights=[2.0, 1.0], low_memory)
        @test R2.md_count ≈ R1.md_count
        @test all(R2.solute_group_count == R1.solute_group_count)
        @test all(R2.solvent_group_count == R1.solvent_group_count)
        R2 = mddf("$data_dir/NAMD/trajectory.dcd", tmao, Options(lastframe=2); frame_weights=[0.5, 0.25], low_memory)
        @test R2.md_count ≈ R1.md_count
        @test all(R2.solute_group_count == R1.solute_group_count)
        @test all(R2.solvent_group_count == R1.solvent_group_count)
        @test R2.volume.total ≈ R1.volume.total
        # Varying weights with stride
        R1 = mddf("$data_dir/NAMD/traj_duplicated_first_frame.dcd", tmao, Options(firstframe=2, lastframe=3); low_memory)
        R2 = mddf("$data_dir/NAMD/traj_duplicated_first_frame.dcd", tmao, Options(lastframe=3, stride=2); frame_weights=[1.0, 100.0, 1.0], low_memory)
        @test R2.md_count ≈ R1.md_count
        @test all(R2.solute_group_count == R1.solute_group_count)
        @test all(R2.solvent_group_count == R1.solvent_group_count)
        @test R2.volume.total ≈ R1.volume.total
    end

end
