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
const bars = "-------------------------------------------------------------------------------"
atoms_str(n) = "$n $(n == 1 ? "atom" : "atoms")"
mol_str(n) = "$n $(n == 1 ? "molecule" : "molecules")"

@testitem "str numbers" begin
    import ComplexMixtures as CM
    @test CM.atoms_str(1) == "1 atom"
    @test CM.atoms_str(2) == "2 atoms"
    @test CM.mol_str(1) == "1 molecule"
    @test CM.mol_str(2) == "2 molecules"
end

#=
    writexyz(x::Vector{T}, file::String) where T <: AbstractVector

Print test xyz file.

=#
function writexyz(x::Vector{T}, file::String) where {T<:AbstractVector}
    f = open(file, "w")
    nx = length(x)
    println(f, nx)
    println(f, "title")
    for i = 1:nx
        println(f, "H $(x[i][1]) $(x[i][2]) $(x[i][3])")
    end
    close(f)
    return nothing
end
