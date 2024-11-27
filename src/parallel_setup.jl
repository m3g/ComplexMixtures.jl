#=

This function will set the number of chunks and the parameters for the CellListMap
parallelization, given the size of the results data structure and the level of parallelization.

=#
function parallel_setup(options::Options, R::Result, low_memory::Bool)

    # Number of threads available 
    nthreads = options.nthreads == 0 ? Threads.nthreads() : options.nthreads

    # Memory requirements, given the size of the results data structure and 
    # the level of parallelization
    total_memory = Sys.total_memory()  / 1024^3
    results_memory = Base.summarysize(R) / 1024^3
    required_memory = results_memory * nthreads

    #
    # Set the number of chunks and CellListMap parameters for low-memory option
    #
    nchunks, parallel_cl, nbatches_cl = if low_memory
        nchunks = max(1, min(nthreads, Int(fld(0.2 * total_memory, results_memory))))
        nthreads_per_chunk = Int(fld(nthreads, nchunks))
        (nchunks, true, (min(8,nthreads_per_chunk), nthreads_per_chunk))
    else
        (nthreads, false, (1,1))
    end

    if nchunks * results_memory > 0.5 * total_memory
        @warn begin 
            """\n
            The memory required for the computation is a large proportion of the total system's memory.
            Depending on resources used by other processes, this may lead to memory exhaustion and
            and the termination of the computation.

            - The Results data structure is $(round(results_memory, digits=2)) GiB, and $nchunks copies are required 
              for the parallel computation, thus requiring $(round(required_memory, digits=2)) GiB.
            - The total system memory is: $(round(total_memory, digits=2)) GiB.

            To reduce memory requirements, consider:

            - Use the `low_memory = true` option in the call to `mddf`.
              For example: `mddf(trajectory_file, solute, solvent, options; low_memory = true)`
            - Using the predefinition of custom groups of atoms in the solute and solvent.
              See: https://m3g.github.io/ComplexMixtures.jl/stable/selection/#predefinition-of-groups
            - Reducing the number of threads, with the `Options(nthreads=N)` parameter.
              Here, we suggest at most N=$(Int(fld(0.2 * total_memory, results_memory))).

            """
        end _file=nothing _line=nothing
    end
    if nchunks * results_memory > total_memory
        throw(ErrorException("The memory required for the computation is larger than the total system memory."))
    end

    return nchunks, parallel_cl, nbatches_cl, nthreads
end
