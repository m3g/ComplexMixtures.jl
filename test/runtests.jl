using TestItemRunner: @run_package_tests, @testitem

using Pkg.Artifacts
testing_dir = artifact"ComplexMixtures"


@run_package_tests

@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(ComplexMixtures)
end
