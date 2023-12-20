using TestItemRunner
@run_package_tests filter = ti -> occursin("AtomSelection.jl", ti.filename)

@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(ComplexMixtures)
end