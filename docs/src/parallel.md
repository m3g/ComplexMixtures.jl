
# Parallel execution 

It is highly recommended to run MDDF calculations in parallel, using multiple processors of a 
single computer. To run the computation in parallel, initialize `julia` with
the `-t N` option, where `N` is the number of processes to be used. For example,
to use `8` parallel processes, use:
```
julia -t 8 example.jl
```
The computation will use a number of parallel processes equal to `N`. Use `-t auto` to automatically
pick the number of threads available in your computer. 

## Optimal number of threads

The number of threads used for computation of the MDDF is the number of threads available to Julia. 
Many computers allow hyperthreading, and not necessarily this this beneficial for the execution
of this package. The optimal number of threads may vary. Some newer CPUs have "energy saving"
cores, which are also relatively slow.

Independently of the number of threads initialized with the `-t` command-line
parameter, the number of processes launched by `ComplexMixtures` in any 
given computation can be adjusted by the `Options(nthreads=N)` option. This
won't provide any speedup if the optional number of threads is greater than
the number of threads available to Julia at runtime.

## Memory issues

If the calculations get `Killed` by no apparent reason, that is probably
because you are running out of memory because of the many parallel computations
running. 

The main reason for memory exhaustion is the annotation of the contributions
of each atom of the solute molecules to the total counts. By default, in
a parallel run, one copy of such data structures is saved for each thread.
This might be an issue if the solute is molecular structure with 
many hundreds of thousand atoms, and if the number of CPUs available is 
very high, but memory is not as abundant.

There are different ways to deal with this issue:

1. Use [predefined groups](@ref predefinition-of-groups), to reduce the
   size of the group array contributions.
2. Use the `low_memory=true` option the call to `mddf` (for example
   with `mddf(traj, options; low_memory=true)`. 
3. Reduce the number of threads used (with `Options(nthreads=N)`).
4. Increase the frequency of garbage collection calls:
   ```julia
   options = Options(GC=true, GC_threshold=0.5)
   R = mddf(trajectory_file, solute, solvent, options)
   ```     
   The `GC_threshold=0.5` indicates that if the free memory is smaller than 50%
   of the total memory of the machine, a garbage-collection run will occur. The  
   default parameters are `GC=true` and `GC_threshold=0.3`.  


