"""
    kbi(R::Result; correction::Symbol=:W7, normalization::Symbol=:ganguly)

Compute the Kirkwood-Buff integral (KBI) from the results of a minimum-distance distribution 
function (MDDF) calculation, `R`, as a function of the upper limit of integration, ``L``. 

Returns a vector with the same length as `R.d`, in cm³ mol⁻¹, where element `i` is the KBI computed 
with ``L`` equal to the upper limit of bin `i`. `R.kb` contains the KBI computed with the default 
options, `R.kb == kbi(R)`.

# Keyword arguments

- `correction`: the weight function that corrects for the truncation of the integral at ``L``:
    - `:W7` (default): the ``W_7^{(3)}`` weight of Santos.
    - `:G2`: the weight of Krüger and Vlugt (Eq. 24).
    - `:G1`: the weight of Krüger et al., ``1 - (d/L)^3``.
    - `:none` (or `:G0`): no weight, the truncated integral.
- `normalization`: the reference density of the solvent:
    - `:ganguly` (default): the density of the solvent outside the domain of the solute within
      each distance, following Ganguly and van der Vegt.
    - `:mddf`: the bulk density of the solvent (estimated beyond the cutoff), at all distances, 
      which is the normalization of `R.mddf`.

The corrections apply to MDDFs, thus to solutes and solvents of any shape, and to radial distribution
functions. Report the value (or average) of the KBI in the range of ``L`` where it is stable. The 
differences between the results obtained with different corrections and normalizations are an 
indication of the uncertainty of the result. 

See the [Kirkwood-Buff integrals](@ref kbi) section of the manual for a description of the 
corrections, references, and examples.

# Examples

```julia-repl
julia> R = load("result.json");

julia> G = kbi(R); # W₇⁽³⁾ weight, Ganguly normalization

julia> G == R.kb
true

julia> G0 = kbi(R; correction=:none, normalization=:mddf); # truncated integral of the MDDF
```

See also [`finite_volume_kbi`](@ref), [`extrapolate_kbi`](@ref), and [`reference_density`](@ref).

!!! compat
    This function was introduced in version 2.19.0. The `:W7` correction and the `normalization` 
    keyword were introduced in version 3.0.0, in which the corrections were extended to MDDFs, 
    and the defaults were changed from `correction=:none` (and the normalization of the MDDF) to 
    `correction=:W7, normalization=:ganguly`. Since version 3.0.0, `R.kb` is equal to `kbi(R)`.

"""
function kbi(R::Result; correction::Symbol=:W7, normalization::Symbol=:ganguly)
    return _kbi(R, R.md_count, R.md_count_random; correction, normalization)
end

#=
    _kbi(R::Result, count, count_random; correction, normalization)

Computes the KBI from the minimum-distance counts `count` and `count_random`, which can 
be those of the complete solute and solvent, or of a group of atoms (in which case the result
is the contribution of the group to the KBI). The reference density is always computed
from the total counts.

=#
function _kbi(R::Result, count::AbstractVector, count_random::AbstractVector; correction::Symbol, normalization::Symbol)
    W = _kbi_weight(Val(correction))
    h = _kbi_integrand(R, count, count_random, Val(normalization))
    return _weighted_running_integral(h, R, W)
end

const _kbi_corrections = (:W7, :G2, :G1, :none, :G0)
function _kbi_weight(::Val{correction}) where {correction}
    throw(ArgumentError("""\n
        Invalid KBI correction option: :$correction. Available corrections are: 
        $(join(":" .* string.(_kbi_corrections), ", "))

    """))
end
# Weights are polynomials, represented by their coefficients (constant term first)
_polymul(a, b) = [sum(a[i] * b[k-i+1] for i in max(1, k - length(b) + 1):min(k, length(a))) for k in 1:(length(a)+length(b)-1)]
_kbi_weight(::Val{:none}) = [1.0]
_kbi_weight(::Val{:G0}) = _kbi_weight(Val(:none))
_kbi_weight(::Val{:G1}) = [1.0, 0.0, 0.0, -1.0]
_kbi_weight(::Val{:G2}) = [1.0, 0.0, 0.0, -23 / 8, 3 / 4, 9 / 8]
# W₇⁽³⁾(x) = (1-x)^4 (1 + 35x/16)(1 + 1225x^2/256)(1 + 29x/16 + 5x^2/4 + 5x^3/16)
_kbi_weight(::Val{:W7}) = foldl(_polymul, (
    [1.0, -4.0, 6.0, -4.0, 1.0], [1.0, 35 / 16], [1.0, 0.0, 1225 / 256], [1.0, 29 / 16, 5 / 4, 5 / 16]
))

const _kbi_normalizations = (:ganguly, :mddf)
function _kbi_integrand(::Result, ::AbstractVector, ::AbstractVector, ::Val{normalization}) where {normalization}
    throw(ArgumentError("""\n
        Invalid KBI normalization option: :$normalization. Available normalizations are: 
        $(join(":" .* string.(_kbi_normalizations), ", "))

    """))
end

function _check_normalized(R::Result)
    if R.density.solvent_bulk == 0
        throw(ArgumentError("""\n
            The bulk solvent density is zero. The Result object does not contain normalized
            distributions (was it computed with `coordination_number`?).

        """))
    end
end

#=
    _kbi_integrand(R::Result, count, count_random, normalization)

Returns, for each bin, the excess volume (in Å³) of the integrand of the KBI, that is, 
the integral of [ρ(d)/ρ_ref(d) - ρ_id(d)/ρ_bulk] dV(d) over the bin, computed from
the minimum-distance counts provided.

=#
function _kbi_integrand(R::Result, count::AbstractVector, count_random::AbstractVector, ::Val{:mddf})
    _check_normalized(R)
    return (count .- count_random) ./ R.density.solvent_bulk
end

function _kbi_integrand(R::Result, count::AbstractVector, count_random::AbstractVector, ::Val{:ganguly})
    _check_normalized(R)
    ρ_ref = _ganguly_density(R)
    return @. count / ρ_ref - count_random / R.density.solvent_bulk
end

#=
    _ganguly_density(R::Result)

Returns the density of the solvent molecules outside the domain of the solute within each
distance `R.d`, ρ_ref(d) = [N - N_in(d)] / [V - V(d)], in molecules/Å³.

=#
function _ganguly_density(R::Result)
    _check_normalized(R)
    Vd = cumsum(R.md_count_random) ./ R.density.solvent_bulk # Volume of the domain within d
    Nin = cumsum(R.md_count) # Number of solvent molecules within d
    N = R.solvent.nmols - (R.autocorrelation ? 1 : 0)
    return @. (N - Nin) / (R.volume.total - Vd)
end

#=
    _weighted_running_integral(h, R::Result, W)

Returns the running integrals ∑ᵢ h[i] W(dᵢ/L) for each upper limit L of the bins, in cm³ mol⁻¹,
where W is a polynomial with coefficients `W` (constant term first). The integrals are computed
from the running moments Sₖ(L) = ∑ᵢ h[i] dᵢᵏ, as ∑ₖ W[k+1] Sₖ(L)/Lᵏ.

=#
function _weighted_running_integral(h::AbstractVector, R::Result, W::AbstractVector)
    binstep = R.files[1].options.binstep
    G = zeros(length(h))
    for k in eachindex(W)
        W[k] == 0 && continue
        Sk = 0.0
        for j in eachindex(G, h)
            Sk += h[j] * R.d[j]^(k - 1)
            G[j] += W[k] * Sk / (j * binstep)^(k - 1)
        end
    end
    return units.Angs3tocm3permol .* G
end

#=
    _rdf_excess(R::Result; caller)

Returns, for each bin, the integral of 4πr²h(r) over the bin, in Å³, computed from the 
RDF counts and the random reference counts. Checks if the RDF-based estimators apply.

=#
function _rdf_excess(R::Result; caller::String)
    _check_normalized(R)
    if R.solute.natomspermol != 1
        @warn begin
            """\n
            `$caller` applies to radial distribution functions, such that the solute 
            selection must contain only one atom per molecule. 

            Here, the number of atoms per molecule of the solute is $(R.solute.natomspermol), implying
            a minimum-distance distribution function, for which the finite-volume interpretation 
            of these estimators (and the 1/L dependence of the KBIs on L) is not justified.

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

The finite-volume theory concerns two-point integrals over open spherical subvolumes, thus ``G(L)`` 
is the KBI of a subvolume, and ``F_\\infty`` a surface term, only for radial distribution functions 
(the solute must contain a single atom per molecule). For an MDDF there is no corresponding spherical 
subvolume: the weighted integral still tends to ``G_\\infty`` as ``L \\to \\infty``, as any weighted 
integral with ``W(0) = 1``, but the ``1/L`` dependence is not justified, and the slope of the 
extrapolation has no physical meaning. For MDDFs, use [`kbi`](@ref).

The weighted estimates of [`kbi`](@ref) do not require extrapolation: the `:G1` weight, 
``W(x) = 1 - x^3``, is ``d[L\\,G(L)]/dL``, which removes the ``1/L`` term of ``G(L)`` analytically.
The extrapolation is useful when ``G(L)`` is available but the running integrals around a 
molecule are not, or as a cross-check of the weighted estimates.

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
    @test_throws "Invalid KBI normalization option" kbi(R; normalization=:abc)
    @test kbi(R; correction=:none, normalization=:mddf) ≈ R.kb
    @test kbi(R; correction=:G0, normalization=:mddf) == kbi(R; correction=:none, normalization=:mddf)
    @test kbi(R) == kbi(R; correction=:W7, normalization=:ganguly)

    # Compare with direct evaluation of the running integrals
    dr = R.files[1].options.binstep
    ρb = R.density.solvent_bulk
    Vd = cumsum(R.md_count_random) ./ ρb
    ρg = (R.solvent.nmols - 1 .- cumsum(R.md_count)) ./ (R.volume.total .- Vd)
    hb = (R.md_count .- R.md_count_random) ./ ρb
    hg = R.md_count ./ ρg .- R.md_count_random ./ ρb
    direct(h, j, w) = units.Angs3tocm3permol * sum(h[i] * w(R.d[i] / (j * dr)) for i in 1:j)
    W7(x) = (1 - x)^4 * (1 + 35x / 16) * (1 + 1225x^2 / 256) * (1 + 29x / 16 + 5x^2 / 4 + 5x^3 / 16)
    g1 = kbi(R; correction=:G1, normalization=:mddf)
    g2 = kbi(R; correction=:G2, normalization=:mddf)
    w7 = kbi(R; correction=:W7, normalization=:mddf)
    w7g = kbi(R)
    fv = finite_volume_kbi(R)
    for j in (100, 500, 1000, length(R.d))
        @test g1[j] ≈ direct(hb, j, x -> 1 - x^3)
        @test g2[j] ≈ direct(hb, j, x -> 1 - 23 / 8 * x^3 + 3 / 4 * x^4 + 9 / 8 * x^5)
        @test w7[j] ≈ direct(hb, j, W7)
        @test w7g[j] ≈ direct(hg, j, W7)
        @test fv.G[j] ≈ direct(hb, j, x -> 1 - 3 / 2 * x + 1 / 2 * x^3)
        @test fv.L[j] ≈ j * dr
    end
    # W7 is the product of y7(x) = (1-x)^4(1 + 29x/16 + 5x^2/4 + 5x^3/16) and ∑ₙ₌₀³ (35x/16)ⁿ
    for x in (0.0, 0.25, 0.5, 0.75, 1.0)
        a = 35 / 16
        @test W7(x) ≈ (1 - x)^4 * (1 + 29x / 16 + 5x^2 / 4 + 5x^3 / 16) * sum((a * x)^n for n in 0:3)
    end

    # Water O-O: G∞ ≈ -16 cm³/mol
    i = findfirst(>=(12.0), R.d):findfirst(>=(20.0), R.d)
    @test all(x -> -17.0 < x < -15.5, g2[i])
    @test all(x -> -17.0 < x < -15.5, w7[i])
    @test all(x -> -17.0 < x < -15.5, w7g[i])
    ext = extrapolate_kbi(fv, (10.0, 20.0))
    @test ext.Ginf ≈ -16.3 atol = 0.2
    @test ext.F ≈ 22.8 atol = 1.0
    @test extrapolate_kbi(R, (10.0, 20.0)) == ext
    # The G1 weight is d[L G(L)]/dL of the finite-volume KBI (removes the 1/L term)
    LG = fv.L .* fv.G
    for j in (300, 500, 1000)
        @test g1[j] ≈ (LG[j+1] - LG[j-1]) / (fv.L[j+1] - fv.L[j-1]) atol = 1e-2
    end

    # Extrapolation of exact data
    L = collect(1.0:0.5:30.0)
    ext = extrapolate_kbi(ComplexMixtures.FiniteVolumeKBI(L, @. 2.0 - 3.0 / L), (5.0, 25.0))
    @test ext.Ginf ≈ 2.0
    @test ext.F ≈ -3.0
    @test_throws ArgumentError extrapolate_kbi(fv, (20.0, 10.0))
    @test_throws ArgumentError extrapolate_kbi(fv, (30.0, 40.0))

    # Multi-atom solute (MDDF): no warning for kbi, warning for finite_volume_kbi
    R = load("$data_dir/NAMD/water/rw_20_25.json")
    @test_logs kbi(R; correction=:G2)
    @test_logs kbi(R)
    @test_logs (:warn,) finite_volume_kbi(R)
    w7g = kbi(R)
    i = findfirst(>=(8.0), R.d):findfirst(>=(12.0), R.d)
    @test all(x -> -17.5 < x < -15.5, w7g[i])

    # Ganguly reference density
    ρg = ComplexMixtures._ganguly_density(R)
    @test R.autocorrelation
    @test ρg[end] ≈ (R.solvent.nmols - 1 - sum(R.md_count)) / (R.volume.total - sum(R.md_count_random) / R.density.solvent_bulk)
    @test all(x -> isapprox(x, R.density.solvent_bulk; rtol=1e-3), ρg[findfirst(>=(5.0), R.d):end])

    # For results computed with the current version, the reference density of the Ganguly 
    # normalization at the cutoff is the bulk density
    R = load("$data_dir/NAMD/protein_tmao.json")
    @test ComplexMixtures._ganguly_density(R)[end] ≈ R.density.solvent_bulk
    R = load("$data_dir/NAMD/tmao_tmao.json")
    @test R.autocorrelation
    @test ComplexMixtures._ganguly_density(R)[end] ≈ R.density.solvent_bulk

    # Not normalized results
    R = load("$data_dir/NAMD/water/rw_20_25.json")
    R.density.solvent_bulk = 0.0
    @test_throws "bulk solvent density is zero" kbi(R)
end
