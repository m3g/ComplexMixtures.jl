#
# Module that provides some global functions for testing
#
module Testing

    const src_dir = @__DIR__
    const test_dir = normpath("$src_dir/../test")
    const data_dir = normpath("$test_dir/data")
    const pdbfile = normpath("$data_dir/NAMD/structure.pdb")

end # module Testing