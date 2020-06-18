#
# Structure that contains the detailed input options
#

@with_kw struct Options

   output :: String

   firstframe :: Int64 = 1
   lastframe :: Int64 = -1
   stride :: Int64 = 1

   periodic :: Bool = true

   irefatom :: Int64 = -1
   n_random_samples :: Int64 = 10

   print_files :: Bool = true 
   print_results :: Bool = true

   binstep :: Float64 = 0.02
   dbulk :: Float64 = 10.
   cutoff :: Float64 = 10.
   usecutoff :: Bool = false

   density_fix :: Bool = false

end

Options() = error(" Options must be initialized at least with the output file name, with Options(output=\"example.dat\")")

Options( output :: String ) = Options(output=output)

