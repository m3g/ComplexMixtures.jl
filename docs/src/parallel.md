
# Parallel execution 

An MDDF calculation can be performed in parallel, using many processors of a 
single computer. The speedup is almost linear, as the parallelization is performed
by splitting the calculation of each frame in each processor, in an
asynchronous manner. To run the computation in parallel, just define the
number of threads (processors) to be used, by defining the
`JULIA_NUM_THREADS` environment variable:
```
export JULIA_NUM_THREADS=4
```
or, which is more simple in Julia 1.5 or greater, just start `julia`
with:  
```
julia -t auto
```

To directly run a script, use
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
