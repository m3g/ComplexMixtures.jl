using TestItemRunner: @run_package_tests, @testitem

@run_package_tests

@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(ComplexMixtures)
end

include("./show_methods.jl")
