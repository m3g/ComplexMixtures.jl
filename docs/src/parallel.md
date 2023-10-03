
# Parallel execution 

It is highly recommended to run MDDF calculations in parallel, using multiple processors of a 
single computer. To run the computation in parallel, initialize `julia` with
the `-t N` option, where `N` is the number of processes to be used. For example,
to use `8` parallel processes, use:
```
julia -t 8 example.jl
```
The computation will use a number of parallel processes equal to `N`. 

!!! note
    The number of threads used for computation of the MDDF is the number of threads available to Julia. 
    Many computers allow hyperthreading, and not necessarily this this beneficial for the execution
    of this package. The optimal number of threads may vary.
    
    Independently of the number of threads initalized with the `-t` command-line
    parameter, the number of processes launched by `ComplexMixtures` in any 
    given computation can be adjusted by the `Options(nthreads=N)` option. This
    won't provide any speedup if the optional number of threads is greater than
    the number of threads available to Julia at runtime.

!!! warning
    If the calculations get `Killed` by no apparent reason, that is probably
    because you are running out of memory because of the many parallel computations
    running. One way to alleviate this problem is to force garbage collection,
    using
    ```julia
    options = Options(GC=true,GC_threshold=0.5)
    R = mddf(trajectory,options)

    ```     
    The `GC_threshold=0.5` indicates that if the free memory is smaller than 50%
    of the total memory of the machine, a garbage-collection run will occur. The  
    default parameters are `GC=true` and `GC_threshold=0.3`.  

