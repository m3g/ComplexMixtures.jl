"""
    kbi(R::Result; correction::Symbol=:none)

Compute the Kirkwood-Buff integral (KBI) from the data in the `Result` object `R`, as a
function of the upper integration limit, optionally using the improved estimators of the
infinite-volume KBI proposed by Krüger and Vlugt [1,2].

The output is a vector with the same length as `R.d`, in units of cm³ mol⁻¹.

# Correction types

- `:none` (default): the KBI computed from the MDDF, without any correction (`R.kb`).
- `:G0`: the KBI computed from the RDF without any correction (`R.kb_rdf`), that is,
  the truncated integral ``G_0(L) = \\int_0^L h(r) 4\\pi r^2 dr``, where ``h(r) = g(r) - 1``.
- `:G1`: the extrapolation of Ref. [1], ``G_1(L) = \\int_0^L h(r) 4\\pi r^2 (1 - x^3) dr``,
  with ``x = r/L``.
- `:G2`: the extrapolation of Ref. [2] (Eq. 24),
  ``G_2(L) = \\int_0^L h(r) 4\\pi r^2 (1 - \\frac{23}{8}x^3 + \\frac{3}{4}x^4 + \\frac{9}{8}x^5) dr``.
  This is the estimator that converges faster to the infinite-volume KBI (with error ``\\sim 1/L^3``).

The `:G1` and `:G2` estimators are computed from the RDF counts (`R.rdf_count` and
`R.rdf_count_random`), and are only valid if the distribution is a radial distribution function,
which requires that the solute contains a single atom per molecule (the solvent can contain
any number of atoms, as the distances are computed to its reference atom).

The estimators correct for the truncation of the integral at a finite distance. They do not
correct the systematic error of the distribution functions computed in closed (canonical)
simulation boxes, which may cause a drift of the KBIs at long distances.

See also [`finite_volume_kbi`](@ref) and [`extrapolate_kbi`](@ref).

# References

1. P. Krüger, S. K. Schnell, D. Bedeaux, S. Kjelstrup, T. J. H. Vlugt, J.-M. Simon,
   Kirkwood–Buff Integrals for Finite Volumes. *J. Phys. Chem. Lett.* 4, 235 (2013).
   [DOI: 10.1021/jz301992u](https://doi.org/10.1021/jz301992u)
2. P. Krüger, T. J. H. Vlugt, Size and shape dependence of finite-volume Kirkwood-Buff integrals.
   *Phys. Rev. E* 97, 051301(R) (2018).
   [DOI: 10.1103/PhysRevE.97.051301](https://doi.org/10.1103/PhysRevE.97.051301)

!!! compat
    This function was introduced in version 2.19.0.

"""
kbi(R::Result; correction::Symbol=:none) = _kbi(R, Val(correction))

function _kbi(::Result, ::Val{correction}) where {correction}
    throw(ArgumentError("""\n
        Invalid KBI correction option: :$correction. Available corrections are:
            - :none
            - :G0
            - :G1
            - :G2

    """))
end

_kbi(R::Result, ::Val{:none}) = copy(R.kb)
_kbi(R::Result, ::Val{:G0}) = copy(R.kb_rdf)
function _kbi(R::Result, ::Val{:G1})
    (; S0, S3, L) = _running_moments(R; caller="kbi")
    return @. units.Angs3tocm3permol * (S0 - S3 / L^3)
end

function _kbi(R::Result, ::Val{:G2})
    (; S0, S3, S4, S5, L) = _running_moments(R; caller="kbi")
    return @. units.Angs3tocm3permol * (S0 - (23 / 8) * S3 / L^3 + (3 / 4) * S4 / L^4 + (9 / 8) * S5 / L^5)
end

#=
    _rdf_excess(R::Result; caller)

Returns, for each bin, the integral of 4πr²h(r) over the bin, in Å³, computed from the 
RDF counts and the random reference counts. Checks if the RDF-based estimators apply.

=#
function _rdf_excess(R::Result; caller::String)
    if R.density.solvent_bulk == 0
        throw(ArgumentError("""\n
            The bulk solvent density is zero. The Result object does not contain normalized
            distributions (was it computed with `coordination_number`?).

        """))
    end
    if R.solute.natomspermol != 1
        @warn begin
            """\n
            `$caller` applies to radial distribution functions, such that the solute 
            selection must contain only one atom per molecule. 

            Here, the number of atoms per molecule of the solute is $(R.solute.natomspermol), implying
            a minimum-distance distribution function, for which these estimators are not valid.

            """
        end _file = nothing _line = nothing
    end
    return (R.rdf_count .- R.rdf_count_random) ./ R.density.solvent_bulk
end

#=
    _running_moments(R::Result; caller)

Returns the running moments Sk(L) = ∫₀ᴸ 4πr²h(r) rᵏ dr (in Å³⁺ᵏ), for k = 0, 1, 3, 4, 5, 
and the upper limits L of each bin. The running integrals with polynomial weights
w(r/L) are linear combinations of these moments.

=#
function _running_moments(R::Result; caller::String)
    h = _rdf_excess(R; caller)
    d = R.d
    L = [i * R.files[1].options.binstep for i in eachindex(d)]
    return (
        S0=cumsum(h),
        S1=cumsum(@. h * d),
        S3=cumsum(@. h * d^3),
        S4=cumsum(@. h * d^4),
        S5=cumsum(@. h * d^5),
        L=L,
    )
end

"""
    FiniteVolumeKBI

Structure that contains the finite-volume KBIs of spheres of diameter `L`, as computed
by [`finite_volume_kbi`](@ref).

# Fields

- `L::Vector{Float64}`: diameters of the spheres (Å).
- `G::Vector{Float64}`: finite-volume KBIs (cm³ mol⁻¹).

!!! compat
    This structure was introduced in version 2.19.0.

"""
struct FiniteVolumeKBI
    L::Vector{Float64}
    G::Vector{Float64}
end

function Base.show(io::IO, ::MIME"text/plain", fv::FiniteVolumeKBI)
    print(io, chomp("""
    FiniteVolumeKBI with $(length(fv.L)) points:
        L from $(round(first(fv.L); digits=3)) to $(round(last(fv.L); digits=3)) Å
        G(L) at the largest L: $(round(last(fv.G); digits=3)) cm³ mol⁻¹
    """))
end

"""
    finite_volume_kbi(R::Result)

Compute the finite-volume Kirkwood-Buff integrals ``G(L)`` of spheres of diameter ``L``,
for each ``L`` up to the cutoff of the distribution function, using the exact weight
function of the sphere [1,2]:

```math
G(L) = \\int_0^L h(r) 4\\pi r^2 \\left(1 - \\frac{3}{2}x + \\frac{1}{2}x^3\\right) dr, \\quad x = r/L
```

where ``h(r) = g(r) - 1``. For large ``L``, ``G(L) = G_\\infty + F_\\infty/L + O(1/L^2)``, where
``G_\\infty`` is the infinite-volume KBI and ``F_\\infty`` is a surface term. Thus, ``G_\\infty``
can be obtained by extrapolating ``G(L)`` as a function of ``1/L`` to ``1/L \\to 0``, which can be
done with the [`extrapolate_kbi`](@ref) function.

Returns a [`FiniteVolumeKBI`](@ref) object, with fields `L` (Å) and `G` (cm³ mol⁻¹).

As with the `:G1` and `:G2` corrections of [`kbi`](@ref), this is only valid for radial
distribution functions, thus the solute must contain a single atom per molecule.

# Example

```julia-repl
julia> R = load("water.json");

julia> fv = finite_volume_kbi(R);

julia> using Plots

julia> plot(1 ./ fv.L, fv.G; xlabel="1/L (Å⁻¹)", ylabel="G(L) (cm³ mol⁻¹)")

julia> extrapolate_kbi(fv, (10.0, 20.0))
(Ginf = -16.3, F = 22.8)
```

# References

1. P. Krüger, S. K. Schnell, D. Bedeaux, S. Kjelstrup, T. J. H. Vlugt, J.-M. Simon,
   Kirkwood–Buff Integrals for Finite Volumes. *J. Phys. Chem. Lett.* 4, 235 (2013).
   [DOI: 10.1021/jz301992u](https://doi.org/10.1021/jz301992u)
2. P. Krüger, T. J. H. Vlugt, Size and shape dependence of finite-volume Kirkwood-Buff integrals.
   *Phys. Rev. E* 97, 051301(R) (2018).
   [DOI: 10.1103/PhysRevE.97.051301](https://doi.org/10.1103/PhysRevE.97.051301)

!!! compat
    This function was introduced in version 2.19.0.

"""
function finite_volume_kbi(R::Result)
    (; S0, S1, S3, L) = _running_moments(R; caller="finite_volume_kbi")
    G = @. units.Angs3tocm3permol * (S0 - (3 / 2) * S1 / L + (1 / 2) * S3 / L^3)
    return FiniteVolumeKBI(L, G)
end

"""
    extrapolate_kbi(fv::FiniteVolumeKBI, Lrange::Tuple{Real,Real})
    extrapolate_kbi(R::Result, Lrange::Tuple{Real,Real})

Extrapolate the finite-volume KBIs ``G(L)`` to the infinite-volume limit, by a linear
least-squares fit of

```math
G(L) = G_\\infty + F_\\infty \\frac{1}{L}
```

to the data with `Lrange[1] <= L <= Lrange[2]` (in Å).

Returns a named tuple `(Ginf, F)`, with the infinite-volume KBI, `Ginf` (cm³ mol⁻¹), and the
surface term, `F` (cm³ mol⁻¹ Å).

The relation is valid only for ``L`` larger than the correlation length of the fluid, and the result
depends on the range chosen. Plot `fv.G` as a function of `1 ./ fv.L` to choose a range
where the dependence is linear.

If a `Result` object is provided, the finite-volume KBIs are computed with [`finite_volume_kbi`](@ref).

!!! compat
    This function was introduced in version 2.19.0.

"""
function extrapolate_kbi(fv::FiniteVolumeKBI, Lrange::Tuple{Real,Real})
    Lmin, Lmax = Lrange
    if !(Lmin < Lmax)
        throw(ArgumentError("The range must satisfy Lmin < Lmax. Got $Lrange."))
    end
    inds = findall(L -> Lmin <= L <= Lmax, fv.L)
    if length(inds) < 2
        throw(ArgumentError("""\n
            Less than two points in the range $Lrange.
            The available values of L range from $(first(fv.L)) to $(last(fv.L)).

        """))
    end
    x = 1 ./ fv.L[inds]
    y = fv.G[inds]
    A = hcat(ones(length(x)), x)
    Ginf, F = A \ y
    return (Ginf=Ginf, F=F)
end
extrapolate_kbi(R::Result, Lrange::Tuple{Real,Real}) = extrapolate_kbi(finite_volume_kbi(R), Lrange)

@testitem "kbi" begin
    using ComplexMixtures
    using ComplexMixtures: data_dir, units
    R = load("$data_dir/NAMD/water/rwO_20_25.json")

    @test_throws "Invalid KBI correction option" kbi(R; correction=:abc)
    @test kbi(R) == R.kb
    @test kbi(R) !== R.kb
    @test kbi(R; correction=:G0) == R.kb_rdf

    # Compare with direct (quadratic) evaluation of the running integrals
    dr = R.files[1].options.binstep
    h = (R.rdf_count .- R.rdf_count_random) ./ R.density.solvent_bulk
    direct(j, w) = units.Angs3tocm3permol * sum(h[i] * w(R.d[i] / (j * dr)) for i in 1:j)
    g0 = kbi(R; correction=:G0)
    g1 = kbi(R; correction=:G1)
    g2 = kbi(R; correction=:G2)
    fv = finite_volume_kbi(R)
    for j in (100, 500, 1000, length(R.d))
        @test g0[j] ≈ direct(j, x -> 1.0)
        @test g1[j] ≈ direct(j, x -> 1 - x^3)
        @test g2[j] ≈ direct(j, x -> 1 - 23 / 8 * x^3 + 3 / 4 * x^4 + 9 / 8 * x^5)
        @test fv.G[j] ≈ direct(j, x -> 1 - 3 / 2 * x + 1 / 2 * x^3)
        @test fv.L[j] ≈ j * dr
    end

    # Water O-O: G∞ ≈ -16 cm³/mol
    i = findfirst(>=(12.0), R.d):findfirst(>=(20.0), R.d)
    @test all(x -> -17.0 < x < -15.5, g2[i])
    ext = extrapolate_kbi(fv, (10.0, 20.0))
    @test ext.Ginf ≈ -16.3 atol = 0.2
    @test ext.F ≈ 22.8 atol = 1.0
    @test extrapolate_kbi(R, (10.0, 20.0)) == ext

    # Extrapolation of exact data
    L = collect(1.0:0.5:30.0)
    ext = extrapolate_kbi(ComplexMixtures.FiniteVolumeKBI(L, @. 2.0 - 3.0 / L), (5.0, 25.0))
    @test ext.Ginf ≈ 2.0
    @test ext.F ≈ -3.0
    @test_throws ArgumentError extrapolate_kbi(fv, (20.0, 10.0))
    @test_throws ArgumentError extrapolate_kbi(fv, (30.0, 40.0))

    # Multi-atom solute: warning
    R = load("$data_dir/NAMD/water/rw_20_25.json")
    @test_logs (:warn,) kbi(R; correction=:G2)
    @test_logs (:warn,) finite_volume_kbi(R)
end
