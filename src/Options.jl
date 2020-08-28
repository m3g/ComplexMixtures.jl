#
# Structure that contains the detailed input options
#

@with_kw struct Options

  firstframe :: Int64 = 1
  lastframe :: Int64 = -1
  stride :: Int64 = 1

  periodic :: Bool = true

  irefatom :: Int64 = -1
  n_random_samples :: Int64 = 10

  binstep :: Float64 = 0.02
  dbulk :: Float64 = 10.
  cutoff :: Float64 = 10.
  usecutoff :: Bool = false

  # Linked cell length will be (cutoff/lcell), might be tunned for maximum
  # performance
  lcell :: Int64 = 2

  # Sleep time between checks of threads for multple spawn
  sleep = 0.01

  # Force garbage collection in parallel runs to avoid memory overflow, 
  # whenever free memory in the system is smaller than GC_threshold
  GC :: Bool = true
  GC_threshold :: Float64 = 0.1

end


