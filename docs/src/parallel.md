
# Parallel execution 

MDDF calculation can be performed in parallel, using many processors of a 
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


