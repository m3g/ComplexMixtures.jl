"""     

```
mddf(trajectory::Trajectory, options::Options)
```

Function that computes the minimum-distance distribution function, atomic contributions, and 
KB integrals, given the `Trajectory` structure of the simulation and, optionally, parameters
given as a second argument of the `Options` type. This is the main function of the `ComplexMixtures` 
package. 

### Examples

```julia-repl
julia> trajectory = ComplexMixtures.Trajectory("./trajectory.dcd",solute,solvent);

julia> results = ComplexMixtures.mddf(trajectory);

```

or

```julia-repl
julia> options = Options(lastframe=1000);

julia> results = ComplexMixtures.mddf(trajectory,options);

```

"""
mddf(trajectory::Trajectory) = mddf(trajectory,Options())

# Choose among serial or parallel version, and self and non-self versions

# In both the self and non-self cases, the number of random samples is 
# n_random_samples*solvent.nmols. However, in the non-self distribution
# the sampling of solvent distances is proportional to the number of solute
# molecules, thus the md count has to be averaged by the solute.nmols. In the
# case of the self-distribution, we compute n(n-1)/2 distances, and we will
# divide this by n random samples, which is the sampling of the random
# distribution. Therefore, we must weight the self-distance count by dividing
# it by (n-1)/2, so that we have a count proportional to n as well, leading
# to the correct weight relative to the random sample. 

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
    samples = Samples(md=(trajectory.solvent.nmols-1)/2,
                      random=options.n_random_samples)
    if nthreads == 1
      mddf_linkedcells(trajectory,options,samples,RNG,mddf_frame_self!)
    else
      mddf_linkedcells_parallel(trajectory,options,samples,RNG,mddf_frame_self!)
    end
  # If solute and solvent are different subsets of the simulation
  else
    samples = Samples(md=trajectory.solute.nmols,
                      random=options.n_random_samples)
    if nthreads == 1
      mddf_linkedcells(trajectory,options,samples,RNG,mddf_frame!)
    else
      mddf_linkedcells_parallel(trajectory,options,samples,RNG,mddf_frame!)
    end
  end
end


