#=

Unit conversions.

=#
@kwdef struct Units{T}
    mole::T = 6.022140857e23
    Angs3tocm3::T = 1e24
    Angs3toL::T = 1e27
    Angs3tocm3permol::T = mole / Angs3tocm3
    Angs3toLpermol::T = mole / Angs3toL
    SitesperAngs3tomolperL::T = Angs3toL / mole
end
const units = Units{Float64}()

#
# Decoration and title
#
const bars = repeat("-", 80)
atoms_str(n) = "$n $(n == 1 ? "atom" : "atoms")"
mol_str(n) = "$n $(n == 1 ? "molecule" : "molecules")"

@testitem "str numbers" begin
    import ComplexMixtures as CM
    @test CM.atoms_str(1) == "1 atom"
    @test CM.atoms_str(2) == "2 atoms"
    @test CM.mol_str(1) == "1 molecule"
    @test CM.mol_str(2) == "2 molecules"
end

function _loading_version_errors(json_version, current_version)
    # Error if the json file is from a newer version than the current one
    if json_version > current_version
        throw(ArgumentError("""\n
            Trying to load a json result file created with a newer version of ComplexMixtures. 
            This can cause unpredictable errors. 

            Current version of ComplexMixtures: $current_version
            Version used to create the output .json file: $json_version

            Please update ComplexMixtures and try again.

        """))
    end
    # Error if a version smaller than 2.0.0 is found
    if json_version < v"2.0.0"
        throw(ArgumentError("""\n
            Trying to load a json result created with an older, and incompatible, version of ComplexMixtures.

            Current version of ComplexMixtures: $current_version
            Version used to create the output .json file: $json_version

            To load this file, install an older version of ComplexMixtures, with, for example:
            
            julia> import Pkg; Pkg.pkg"add ComplexMixtures@$json_version"

            You can pin the version of ComplexMixtures to the one you installed with:

            julia> import Pkg; Pkg.pkg"pin ComplexMixtures@$json_version"

        """))
    end
    return nothing
end

# function that checks if output files were produced with the current version
function _check_version(filename)
    json_version = _get_version(filename)
    current_version = pkgversion(@__MODULE__)
    _loading_version_errors(json_version, current_version)
    #
    # Choose load type depending on the version of the file
    #
    load_type = if v"2.0.0" < json_version <= v"2.17.1"
        @warn """\n
            Reading json result created with $json_version. 
            
            This version didn't support KBI contribution decomposition, introduced in v2.18.0. 
            
            The `solute_group_count_random` and `solvent_group_count_random` arrays will be empty.

        """ _file=nothing _line=nothing
        Result_2_17_1
    else
        Result
    end
    return load_type
end