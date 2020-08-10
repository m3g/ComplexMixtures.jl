#     
# mddf_choose
#
# Select serial of parallel version depending on the number of available
# threads
#  

# With default input options

mddf_choose(trajectory) = mddf_choose(trajectory,Options())

# Choose among serial or parallel version

function mddf_choose(trajectory,options)
  nthreads = Threads.nthreads()
  if nthreads == 1
    mddf_linkedcells(trajectory,options)
  else
    mddf_linkedcells_parallel(trajectory,options)
  end
end

