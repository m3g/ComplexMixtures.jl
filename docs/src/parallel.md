
# Parallel execution 

It is highly recommended to run MDDF calculations in parallel, using multiple processors of a 
single computer. To run the computation in parallel, initialize `julia` with
the `-t auto` option:
```
julia -t auto
```
The computation will use a number of threads equal to the number
of physical cores of the computer. The number of computation threads to 
be used can be set by the `Options(nthreads=N)` parameter, where `N` is
an integer. Hyperthreading (using more threads than physical CPUs) 
usually does not provide a significant speedup, and can be detrimental 
in some cases.  

To directly run a script in parallel, use:
```julia
julia -t auto example.jl
```

!!! note
    The number of threads used for computation of the MDDF is the number of threads available to Julia. 
    Many computers allow hyperthreading, and not necessarily this this beneficial for the execution
    of this package. The optimal number of threads may vary.
    
    To adjust the number of threads use the `Options(nthreads=N)` parameter.

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
    default parameters are `GC=true` and `GC_threshold=0.1`.  

