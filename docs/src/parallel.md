
# Parallel execution 

An MDDF calculation can be performed in parallel, using many processors of a 
single computer. To run the computation in parallel, initialize `julia` with
the `-t auto` option:
```
julia -t auto
```
The computation will a number of computation threads equal to the number
of physical cores of the computer. The number of computation threads to 
be used can be set by the `Option(nthreads=N)` parameter, where `N` is
an integer. Hyperthreading (using more threads than physical CPUs) 
usually does not provide a significant speedup, and can be detrimental 
in some cases.  

To directly run a script in parallel, use:
```julia
julia -t auto example.jl
```

!!! note
    The number of threads used for computation of the MDDF is the number
    of physical CPUs of the computer, which are obtained programmatically.
    Most times the use of hyper-threading is not beneficial here.

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

    Unfortunately, this may slow the calculations quite a bit, and the parallelization
    to many processors becomes not very satisfactory. We are working to improve
    this.
