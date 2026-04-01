"""
    renormalize(R::Result, bulk_density::Number, unit::String="mol/L"; silent=true)

Renormalizes the `Result` structure `R` to a different bulk density of the solvent.
The `unit` argument can be either `"mol/L"` or `"sites/Angs3"` (default is `"mol/L"`).

This function does not modify the input Result structure, but returns a new one.

# Example

In the following example we multiply the density of the bulk solvent by 2, causing
the decrease of the MDDF by a factor of 2. The KB integrals will be also updated 
accordingly.

```jldoctest
julia> using ComplexMixtures

julia> R = load(ComplexMixtures.data_dir*"/NAMD/protein_tmao.json");

julia> R_new = renormalize(R, 2*R.density.solvent_bulk, "sites/Angs3");

julia> R_new.mddf ≈ 0.5 * R.mddf
true
```

The `silent` keyword argument controls whether warnings are printed associated to
bins with zero samples in the ideal-gas histogram.

"""
function renormalize(R::Result, bulk_density::Number, unit::String="mol/L"; silent=true)
    if !(unit in ("mol/L", "sites/Angs3"))
        throw(ArgumentError("""\n
            Concentration unit must be "mol/L" or "sites/Angs3".

        """))
    end
    if unit == "mol/L"
        bulk_density = bulk_density / units.SitesperAngs3tomolperL  
    end
    R_new = deepcopy(R) # better define a custom copy function
    density_fix = bulk_density / R.density.solvent_bulk 
    return renormalize!(R_new, density_fix; silent)
end

@testitem "renormalize" begin
    using ComplexMixtures 
    using ComplexMixtures: data_dir, units
    R = load(joinpath(data_dir,"NAMD/protein_tmao.json"))
    @test R.density.solvent != R.density.solvent_bulk
    new_density = 2 * R.density.solvent_bulk * units.SitesperAngs3tomolperL
    R_new = renormalize(R, new_density)
    @test R_new.mddf ≈ 0.5 * R.mddf
    @test R_new.rdf ≈ 0.5 * R.rdf
    @test contributions(R, SolventGroup(["H11", "H12", "H13", "H21", "H22", "H23"])) ≈ 
          2*contributions(R_new, SolventGroup(["H11", "H12", "H13", "H21", "H22", "H23"]))
    @test R_new.kb ≈ 
        units.Angs3tocm3permol * (1/R_new.density.solvent_bulk) * 
        (((R.density.solvent_bulk/units.Angs3tocm3permol) * R.kb) .- R.coordination_number_random)
    @test R_new.kb_rdf ≈ 
        units.Angs3tocm3permol * (1/R_new.density.solvent_bulk) * 
        (((R.density.solvent_bulk/units.Angs3tocm3permol) * R.kb_rdf) .- R.sum_rdf_count_random)

    new_density = 2 * R.density.solvent_bulk
    R_new = renormalize(R, new_density, "sites/Angs3")
    @test R_new.mddf ≈ 0.5 * R.mddf
    @test R_new.rdf ≈ 0.5 * R.rdf
    @test R_new.kb ≈ 
        units.Angs3tocm3permol * (1/R_new.density.solvent_bulk) * 
        (((R.density.solvent_bulk/units.Angs3tocm3permol) * R.kb) .- R.coordination_number_random)
    @test R_new.kb_rdf ≈ 
        units.Angs3tocm3permol * (1/R_new.density.solvent_bulk) * 
        (((R.density.solvent_bulk/units.Angs3tocm3permol) * R.kb_rdf) .- R.sum_rdf_count_random)

    @test_throws ArgumentError renormalize(R, new_density, "wrong_units") 
end