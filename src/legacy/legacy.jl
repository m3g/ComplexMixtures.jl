# some legacy interfaces, to not break old code

import Base: getproperty
getproperty(r::Result, s::Symbol) = getproperty(r, Val(s))
getproperty(r::Result, ::Val{S}) where {S} = getfield(r, S)
getproperty(r::Result, ::Val{:sum_md_count}) = r.coordination_number
getproperty(r::Result, ::Val{:sum_md_count_random}) = r.coordination_number_random

export contrib
const contrib = contributions

@testitem "legacy" begin
    using ComplexMixtures, PDBTools
    import ComplexMixtures.Testing: test_dir

    pdb = readPDB("$test_dir/data/NAMD/structure.pdb")
    R = load("$test_dir/data/NAMD/protein_tmao.json")

    @test R.sum_md_count == R.coordination_number
    @test R.sum_md_count_random == R.coordination_number_random

    @test ComplexMixtures.contrib == ComplexMixtures.contributions
end

#
# Function to load legacy result data structures
#
include("./data_structs_1.0.0-1.3.4.jl")

function load_legacy_json(filename, json_version)
    if json_version >= v"1.0.0" && json_version <= v"1.3.4"
        results = load_json_LegacyA(filename)
    end
    return results
end

