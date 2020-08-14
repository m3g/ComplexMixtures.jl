
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
julia -t 4 
```

To directly run a script, use
```julia
julia -t 4 example.jl

```

!!! note
    The number of threads used for computation of the MDDF is the number
    of defined threads minus one, because one thread is dedicated to
    control the execution. Since the control of the execution is not 
    very demanding, particularly if the number of threads is small, 
    you may want to set the number of threads as one 
    more than you originally intended, if the total number of threads
    is not very large. In particular, running with only `-t 2` will 
    not parallelize the calculation at all.

!!! warning
    If the calculations get `Killed` by no apparent reason, that is probably
    because you are running out of memory because of the many parallel computations
    running. One way to aleviate this problem is to force garbage collection,
    using
    ```julia
    options = ComplexMixtures.Options(GC=true)
    R = ComplexMixtures.mddf(trajectory,options)

    ```     
    Unfortunately, this slows the calculations quite a bit, and the parallelization
    to many processors becomes not very satisfactory. We are working to improve
    this.


