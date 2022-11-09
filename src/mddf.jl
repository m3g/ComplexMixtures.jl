"""     
    mddf(trajectory::Trajectory, options::Options)

Function that computes the minimum-distance distribution function, atomic contributions, and 
KB integrals, given the `Trajectory` structure of the simulation and, optionally, parameters
given as a second argument of the `Options` type. This is the main function of the `ComplexMixtures` 
package. 

### Examples

```julia-repl
julia> trajectory = Trajectory("./trajectory.dcd",solute,solvent);

julia> results = mddf(trajectory);
```

or, to set some custom optional parameter,

```julia-repl
julia> options = Options(lastframe=1000);

julia> results = mddf(trajectory,options);
```

"""
mddf(trajectory::Trajectory) = mddf(trajectory, Options())


function mddf(trajectory::Trajectory, options::Options)
    # Set random number generator
    RNG = init_random(options)
    if options.nthreads < 0
        nthreads = Threads.nthreads()
    else
        nthreads = options.nthreads
    end
    # If the solute and the solvent are the same
    if trajectory.solute.index == trajectory.solvent.index
        samples = Samples(md = (trajectory.solvent.nmols - 1) / 2, random = options.n_random_samples)
        comp_type = Self()
    # If solute and solvent are different subsets of the simulation
    else
        samples = Samples(md = trajectory.solute.nmols, random = options.n_random_samples)
        comp_type = Cross()
    end
    R = mddf_compute(comp_type,trajectory, options, samples, RNG)
    return R
end

"""     
    mddf_compute(trajectory::Trajectory, options::Options, samples::Samples, RNG)  

$(INTERNAL)

Computes the MDDFs. Returns the `Result` structure. This function is resposible for dispatching
the computations of each frame to each thread.

"""
function mddf_compute(comp_type,trajectory::Trajectory, options::Options, samples::Samples, RNG)

    # Number of threads (chunks) to use
    nchunks = options.nthreads
    if nchunks <= 0
        nchunks = Threads.nthreads()
    end

    # Structure in which the final results will be stored
    R0 = Result(trajectory, options)

    # Initializing the structure that carries the result per thread
    R = [Result(trajectory, options) for _ in 1:nchunks]

    # Initialize periodic system for CellListMap
    system = setup_PeriodicSystem(comp_type, trajectory, options)    

    # Create data structures required for multithreading: needed to read coordinates in each
    # frame independently, and compute the minimum-distance list 
    system_chunk = [ deepcopy(system) for _ in 1:nchunks ]

    # Skip initial frames if desired
    iframe = 0 # Counter for all frames of the input file, that must be read serially
    while iframe < options.firstframe - 1
        nextframe!(trajectory)
        iframe += 1
    end

    # Print some information about this run
    options.silent || title(R[1], trajectory.solute, trajectory.solvent, nchunks)

    if !options.silent
        progress = Progress(R0.nframes_read, 1)
    end
    read_lock = ReentrantLock()
    Threads.@threads for (frame_range, ichunk) in chunks(1:R0.nframes_read, nchunks)
        for _ in frame_range
            # Read frame coordinates
            lock(read_lock) do 
                iframe += 1
                while (iframe+1)%options.stride != 0
                    nextframe!(trajectory)
                    iframe += 1
                end
                @. system_chunk[ichunk].x_positions .= trajectory.x_solute
                @. system_chunk[ichunk].y_positions .= trajectory.x_solvent
                update_unitcell!(system_chunk[ichunk], getsides(trajectory, iframe))
                # Run GC if memory is getting full: this are issues with Chemfiles reading scheme
                if options.GC && (Sys.free_memory() / Sys.total_memory() < options.GC_threshold)
                    GC.gc() 
                end
                next!(progress)
            end # reading lock
            # Compute distances in this frame and update results
            mddf_frame!(comp_type, system_chunk[ichunk], options, RNG, R[ichunk])
        end
    end
    closetraj(trajectory)

    # Sum up the results of all threads into the data of thread one (R1<-R1+R2)
    for ichunk in 1:nchunks
        R0 = sum!(R0, R[ichunk])
    end

    # Setup the final data structure with final values averaged over the number of frames,
    # sampling, etc, and computes final distributions and integrals
    R0 = finalresults!(R0, options, trajectory, samples)
    options.silent || println(bars)

    return R0
end

"""
    mddf_frame(::Cross, iframe::Int, framedata::FrameData, options::Options, RNG, R::Result, system::PeriodicSystem)

$(INTERNAL)

Computes the MDDF for a single frame, modifies the data in the `R` (type `Result`) structure.

"""
function mddf_frame!(::Cross, iframe::Int, framedata::FrameData, options::Options, RNG, R::Result, system::PeriodicSystem)

    # Simplify code by assigning some shortened names
    trajectory = framedata.trajectory
    volume_frame = framedata.volume_frame
    rdf_count_random_frame = framedata.rdf_count_random_frame
    md_count_random_frame = framedata.md_count_random_frame
    dc = framedata.dc
    dmin_mol = framedata.dmin_mol
    dref_mol = framedata.dref_mol
    x_solvent_random = framedata.x_solvent_random
    lc_solvent = framedata.lc_solvent
    solute = trajectory.solute
    solvent = trajectory.solvent
    x_solute = trajectory.x_solute
    x_solvent = trajectory.x_solvent

    # Reset counters for this frame
    reset!(volume_frame)
    @. rdf_count_random_frame = 0.0
    @. md_count_random_frame = 0.0

    # Check if the cutoff is not too large considering the periodic cell size
    if R.cutoff > sides[1] / 2 || R.cutoff > sides[2] / 2 || R.cutoff > sides[3] / 2
        error("""
        in MDDF: cutoff or dbulk > periodic_dimension/2 in frame: $iframe
                 max(cutoff,dbulk) = $(R.cutoff)
                 sides read from file = $(sides)
                 If sides are zero it is likely that the PBC information is not available 
                 in the trajectory file.
        """)
    end

    volume_frame.total = sides[1] * sides[2] * sides[3]
    R.volume.total = R.volume.total + volume_frame.total

    R.density.solute = R.density.solute + (solute.nmols / volume_frame.total)
    R.density.solvent = R.density.solvent + (solvent.nmols / volume_frame.total)

    # Compute minimum distances of the molecules to the solute
    minimum_distances!(system)

    # Fraction of solvent molecules in bulk
    n_solvent_in_bulk = count(mol -> !mol.within_cutoff, list) / solute.nmols

    local n_dmin_in_bulk
    n_solvent_in_bulk = 0.0
    for isolute = 1:solute.nmols

        # We need to do this one solute molecule at a time to avoid exploding the memory
        # requirements
        x_this_solute = viewmol(isolute, x_solute, solute)

        # Compute all distances between solute and solvent atoms which are smaller than the 
        # cutoff (this is the most computationally expensive part), the distances are returned
        # in the dc structure
        cutoffdistances!(R.cutoff, x_this_solute, x_solvent, lc_solvent, box, dc)

        # For each solute molecule, update the counters (this is highly suboptimal, because
        # within updatecounters there are loops over solvent molecules, in such a way that
        # this will loop with cost nsolute*nsolvent. However, I cannot see an easy solution 
        # at this point with acceptable memory requirements
        n_dmin_in_bulk, n_dref_in_bulk = updatecounters!(R, solute, solvent, dc, dmin_mol, dref_mol)
        n_solvent_in_bulk += n_dref_in_bulk
    end
    n_solvent_in_bulk = n_solvent_in_bulk / solute.nmols

    #
    # Computing the random-solvent distribution to compute the random minimum-distance count
    #
    bulk_range = (solvent.nmols-n_dmin_in_bulk+1):solvent.nmols
    for i = 1:options.n_random_samples

        # generate random solvent box, and store it in x_solvent_random
        for j = 1:solvent.nmols
            # Choose randomly one molecule from the bulk, if there are actual bulk molecules
            if n_dmin_in_bulk != 0
                jmol = dmin_mol[random(RNG, bulk_range)].jmol
            else
                jmol = random(RNG, 1:solvent.nmols)
            end
            # Indexes of this molecule in the x_solvent array
            x_ref = viewmol(jmol, x_solvent, solvent)
            # Indexes of the random molecule in random array
            x_rnd = viewmol(j, x_solvent_random, solvent)
            # Generate new random coordinates (translation and rotation) for this molecule
            random_move!(x_ref, R.irefatom, sides, x_rnd, RNG)
        end

        # Choose randomly one solute molecule to be the solute in this sample
        i_rand_mol = random(RNG, 1:solute.nmols)
        x_this_solute = viewmol(i_rand_mol, x_solute, solute)

        # Compute all distances between solute and solvent atoms which are smaller than the 
        # cutoff (this is the most computationally expensive part), the distances are returned
        # in the dc structure
        cutoffdistances!(R.cutoff, x_this_solute, x_solvent_random, lc_solvent, box, dc)

        # Update the counters of the random distribution
        updatecounters!(R, rdf_count_random_frame, md_count_random_frame, solvent, dc, dmin_mol, dref_mol)

    end # random solvent sampling

    # Update counters with the data of this frame
    update_counters_frame!(R, rdf_count_random_frame, md_count_random_frame, volume_frame, solute, solvent, n_solvent_in_bulk)

    return nothing
end

"""
    mddf_frame!(::Self, iframe::Int, framedata::FrameData, options::Options, RNG, R::Result)

$(INTERNAL)

Computes the MDDF for a single frame, modifies the data in the `R` (`Result`) structure.

"""
function mddf_frame!(::Self, iframe::Int, framedata::FrameData, options::Options, RNG, R::Result)

    # Simplify code by assigning some shortened names
    trajectory = framedata.trajectory
    volume_frame = framedata.volume_frame
    rdf_count_random_frame = framedata.rdf_count_random_frame
    md_count_random_frame = framedata.md_count_random_frame
    dc = framedata.dc
    dmin_mol = framedata.dmin_mol
    dref_mol = framedata.dref_mol
    x_solvent_random = framedata.x_solvent_random
    lc_solvent = framedata.lc_solvent
    solute = trajectory.solute
    solvent = trajectory.solvent
    x_solute = trajectory.x_solute
    x_solvent = trajectory.x_solvent

    # Reset counters for this frame
    reset!(volume_frame)
    @. rdf_count_random_frame = 0.0
    @. md_count_random_frame = 0.0

    # get pbc sides in this frame
    sides = getsides(trajectory, iframe)

    volume_frame.total = sides[1] * sides[2] * sides[3]
    R.volume.total = R.volume.total + volume_frame.total

    R.density.solute = R.density.solute + (solute.nmols / volume_frame.total)
    R.density.solvent = R.density.solvent + (solvent.nmols / volume_frame.total)

    # Add the box side information to the box structure, in this frame
    box = Box(options.lcell, sides, R.cutoff)

    # Will wrap everything relative to the reference atom of the first molecule
    # and move everything such that that center is in the origin. This is important
    # to simplify the computation of cell indexes, as the minimum coordinates are 
    # automatically -side/2 at each direction
    solute_center = x_solute[R.irefatom]
    wrap!(x_solvent, sides, solute_center)
    center_to_origin!(x_solvent, solute_center)

    # Initialize linked cells
    initcells!(x_solvent, box, lc_solvent)

    # Check if the cutoff is not too large considering the periodic cell size
    if R.cutoff > sides[1] / 2.0 || R.cutoff > sides[2] / 2.0 || R.cutoff > sides[3] / 2.0
        error(
            """
        in MDDF: cutoff or dbulk > periodic_dimension/2 in frame: $iframe
                 max(cutoff,dbulk) = $(R.cutoff)
                 sides read from file = $(sides)
                 If sides are zero it is likely that the PBC information is not available 
                 in the trajectory file.
      """,
        )
    end

    n_solvent_in_bulk = 0.0
    local n_dmin_in_bulk
    for isolvent = 1:solvent.nmols-1
        # We need to do this one solute molecule at a time to avoid exploding the memory
        # requirements
        x_this_solute = viewmol(isolvent, x_solvent, solvent)

        # Compute all distances between solute and solvent atoms which are smaller than the 
        # cutoff (this is the most computationally expensive part), the distances are returned
        # in the dc structure
        cutoffdistances_self!(
            R.cutoff,
            x_this_solute,
            x_solvent,
            lc_solvent,
            box,
            dc,
            solvent,
            isolvent,
        )

        # For each solute molecule, update the counters (this is highly suboptimal, because
        # within updatecounters there are loops over solvent molecules, in such a way that
        # this will loop with cost nsolute*nsolvent. However, I cannot see an easy solution 
        # at this point with acceptable memory requirements
        n_dmin_in_bulk, n_dref_in_bulk =
            updatecounters!(R, solute, solvent, dc, dmin_mol, dref_mol)
        n_solvent_in_bulk += n_dref_in_bulk
    end
    # The normalization below is tricky. The number that comes out from updatecounters is the
    # sum, for every solvent molecule (minus one) of the distances that were not found to be in
    # the solute domain. The total number of distances is n^2, because the sum inside updatecounter
    # is made for all molecules (the same function is used for cross-distribution), but we called updatecounters
    # only (n-1) times, so the actual sum of the distances considered is n(n-1). From this set
    # the number of distances returned in n_dref_in_bulk is r=n(n-1)-nd/2, where nd is the
    # number of distances in the domain. Thus, the number of distances in the domain, considering
    # symmetric terms, is nd=2n(n-1)-2r. The average number of distances in the domain, per
    # molecule, is thus nd/n=2(n-1)-2r/n. Finally, the number of solvent molecules in 
    # bulk, for each molecule, is the total number of other molecules, (n-1), minus the
    # number of molecules in the domain, that is (n-1)-nd/n=(n-1)-2(n-1)+2r/n, which
    # finally simplifies to 2r/n-(n-1), which is the equation below. 
    n_solvent_in_bulk = 2 * n_solvent_in_bulk / solvent.nmols - (solvent.nmols - 1)

    #
    # Computing the random-solvent distribution to compute the random minimum-distance count
    #
    for i = 1:options.n_random_samples
        # generate random solvent box, and store it in x_solvent_random
        for j = 1:solvent.nmols
            # Choose randomly one molecule from the bulk, if there are actually bulk molecules
            if n_dmin_in_bulk != 0
                jmol =
                    dmin_mol[random(RNG, solvent.nmols-n_dmin_in_bulk+1:solvent.nmols)].jmol
            else
                jmol = random(RNG, 1:solvent.nmols)
            end
            # Indexes of this molecule in the x_solvent array
            x_ref = viewmol(jmol, x_solvent, solvent)
            # Indexes of the random molecule in random array
            x_rnd = viewmol(j, x_solvent_random, solvent)
            # Generate new random coordinates (translation and rotation) for this molecule
            random_move!(x_ref, R.irefatom, sides, x_rnd, RNG)
        end

        # wrap random solvent coordinates to origin
        wrap!(x_solvent_random, sides)

        # Initialize linked cells
        initcells!(x_solvent_random, box, lc_solvent)

        # Choose randomly one solute molecule to be the solute in this sample
        i_rand_mol = random(RNG, 1:solvent.nmols)
        x_this_solute = viewmol(i_rand_mol, x_solvent, solvent)

        # Compute all distances between solute and solvent atoms which are smaller than the 
        # cutoff (this is the most computationally expensive part), the distances are returned
        # in the dc structure (here we do not use the self function)
        cutoffdistances!(R.cutoff, x_this_solute, x_solvent_random, lc_solvent, box, dc)

        # Update the counters and get the number of solvent molecules in bulk
        updatecounters!(
            R,
            rdf_count_random_frame,
            md_count_random_frame,
            solvent,
            dc,
            dmin_mol,
            dref_mol,
        )

    end # random solvent sampling

    # Update global counters with the data of this frame
    update_counters_frame!(
        R,
        rdf_count_random_frame,
        md_count_random_frame,
        volume_frame,
        solute,
        solvent,
        n_solvent_in_bulk,
    )

    return R
end







