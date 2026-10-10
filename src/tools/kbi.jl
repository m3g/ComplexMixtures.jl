"""
    kbi(R::Result; correction::Symbol=:W7, normalization::Symbol=:ganguly)

Compute the Kirkwood-Buff integral (KBI) from the minimum-distance distribution function (MDDF)
in the `Result` object `R`, as a function of the upper integration limit ``L``. This is the 
recommended way to obtain the KBI from a `Result`.

The output is a vector with the same length as `R.d`, in units of cm³ mol⁻¹, where the element 
`i` is the KBI computed with ``L`` equal to the upper limit of bin `i`.

The KBI is computed as

```math
G(L) = \\int_0^L \\left[\\frac{\\rho(d)}{\\rho_{\\rm ref}(d)} - \\frac{\\rho_{\\rm id}(d)}{\\rho_{\\rm bulk}}\\right] W(d/L)\\, dV(d)
```

where ``\\rho(d)`` and ``\\rho_{\\rm id}(d)`` are the densities of the solvent at distance ``d`` of the
solute, in the simulation and in the ideal-gas reference, ``W(x)`` is a weight function that corrects 
for the truncation of the integral at a finite distance (``W(x) = 1`` for the truncated integral), 
and ``\\rho_{\\rm ref}(d)`` is the reference density.

# Corrections (weight functions)

- `:W7` (default): the weight ``W_7^{(3)}(x) = (1-x)^4 (1 + \\frac{35}{16}x)(1 + \\frac{1225}{256}x^2)(1 + \\frac{29}{16}x + \\frac{5}{4}x^2 + \\frac{5}{16}x^3)``
  recommended by Santos [3], which vanishes smoothly at ``x = 1``.
- `:G2`: the weight of Krüger and Vlugt [2] (Eq. 24), ``W(x) = 1 - \\frac{23}{8}x^3 + \\frac{3}{4}x^4 + \\frac{9}{8}x^5``.
- `:G1`: the weight of Krüger et al. [1], ``W(x) = 1 - x^3``.
- `:none` (or `:G0`): the truncated integral, ``W(x) = 1``.

The weight functions were derived for radial distribution functions, but, as shown by Santos [3], 
they follow from a purely mathematical identity, valid for any one-dimensional integral. Thus, 
they can be applied to the integrand of the KBI computed from MDDFs. The weights correct the 
truncation of the integral when the integrand oscillates around zero at ``L``. They do not correct 
slowly decaying (monotonic) tails of the distribution, nor errors in the reference density.

# Normalizations (reference density)

- `:ganguly` (default): the reference density at each distance is the density of the solvent 
  molecules outside the domain of the solute within ``d``, as proposed by Ganguly and van der Vegt [4]:
  ``\\rho_{\\rm ref}(d) = [N - N_{\\rm in}(d)]/[V - V(d)]``, where ``N`` is the number of solvent
  molecules (minus one if the solute and solvent are the same), ``N_{\\rm in}(d)`` is the number of 
  solvent molecules within ``d``, and ``V(d)`` is the volume of the domain within ``d``. This corrects 
  for the depletion (or excess) of solvent molecules in the bulk of a closed simulation box.
- `:bulk`: the reference density is the bulk density of the solvent, `R.density.solvent_bulk`, 
  estimated in the bulk region defined by the `dbulk` and `cutoff` parameters. With this normalization 
  and `correction=:none`, the result is `R.kb`.

The normalizations differ by the region of the system from which the reference density is estimated: 
the region beyond ``d`` (`:ganguly`), or the bulk region (`:bulk`). If they lead to different KBIs, the 
estimate of the reference density is a source of error that must be considered.

# Examples

```julia-repl
julia> R = load("result.json");

julia> G = kbi(R); # W₇⁽³⁾ weight, Ganguly normalization

julia> G0 = kbi(R; correction=:none, normalization=:bulk); # equal to R.kb
```

# References

1. P. Krüger, S. K. Schnell, D. Bedeaux, S. Kjelstrup, T. J. H. Vlugt, J.-M. Simon,
   Kirkwood–Buff Integrals for Finite Volumes. *J. Phys. Chem. Lett.* 4, 235 (2013).
   [DOI: 10.1021/jz301992u](https://doi.org/10.1021/jz301992u)
2. P. Krüger, T. J. H. Vlugt, Size and shape dependence of finite-volume Kirkwood-Buff integrals.
   *Phys. Rev. E* 97, 051301(R) (2018).
   [DOI: 10.1103/PhysRevE.97.051301](https://doi.org/10.1103/PhysRevE.97.051301)
3. A. Santos, Finite-size estimates of Kirkwood-Buff and similar integrals.
   [arXiv:1806.00821](https://arxiv.org/abs/1806.00821) (2018).
4. P. Ganguly, N. F. A. van der Vegt, Convergence of Sampling Kirkwood–Buff Integrals of Aqueous Solutions 
   with Molecular Dynamics Simulations. *J. Chem. Theory Comput.* 9, 1347 (2013).
   [DOI: 10.1021/ct301017q](https://doi.org/10.1021/ct301017q)

See also [`finite_volume_kbi`](@ref) and [`extrapolate_kbi`](@ref).

!!! compat
    This function was introduced in version 2.19.0. The `:W7` correction and the `normalization` 
    keyword were introduced in version 2.19.1, in which the corrections were extended to MDDFs, 
    and the defaults were changed from `correction=:none` (and the bulk normalization) to 
    `correction=:W7, normalization=:ganguly`.

"""
function kbi(R::Result; correction::Symbol=:W7, normalization::Symbol=:ganguly)
    W = _kbi_weight(Val(correction))
    h = _kbi_integrand(R, Val(normalization))
    return _weighted_running_integral(h, R, W)
end

const _kbi_corrections = (:W7, :G2, :G1, :none, :G0)
function _kbi_weight(::Val{correction}) where {correction}
    throw(ArgumentError("""\n
        Invalid KBI correction option: :$correction. Available corrections are: 
        $(join(":" .* string.(_kbi_corrections), ", "))

    """))
end
_kbi_weight(::Val{:none}) = x -> 1.0
_kbi_weight(::Val{:G0}) = _kbi_weight(Val(:none))
_kbi_weight(::Val{:G1}) = x -> 1 - x^3
_kbi_weight(::Val{:G2}) = x -> 1 - (23 / 8) * x^3 + (3 / 4) * x^4 + (9 / 8) * x^5
_kbi_weight(::Val{:W7}) = x -> (1 - x)^4 * (1 + 35x / 16) * (1 + 1225x^2 / 256) * (1 + 29x / 16 + 5x^2 / 4 + 5x^3 / 16)

const _kbi_normalizations = (:ganguly, :bulk)
function _kbi_integrand(::Result, ::Val{normalization}) where {normalization}
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
    _kbi_integrand(R::Result, normalization)

Returns, for each bin, the excess volume (in Å³) of the integrand of the KBI, that is, 
the integral of [ρ(d)/ρ_ref(d) - ρ_id(d)/ρ_bulk] dV(d) over the bin.

=#
function _kbi_integrand(R::Result, ::Val{:bulk})
    _check_normalized(R)
    return (R.md_count .- R.md_count_random) ./ R.density.solvent_bulk
end

function _kbi_integrand(R::Result, ::Val{:ganguly})
    _check_normalized(R)
    ρ_ref = _ganguly_density(R)
    return @. R.md_count / ρ_ref - R.md_count_random / R.density.solvent_bulk
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

Returns the running integrals ∑ᵢ h[i] W(dᵢ/L) for each upper limit L of the bins, in cm³ mol⁻¹.

=#
function _weighted_running_integral(h::AbstractVector, R::Result, W::F) where {F<:Function}
    binstep = R.files[1].options.binstep
    G = zeros(length(h))
    for j in eachindex(G)
        L = j * binstep
        G[j] = units.Angs3tocm3permol * sum(h[i] * W(R.d[i] / L) for i in 1:j)
    end
    return G
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

This is only valid for radial distribution functions, thus the solute must contain a single atom 
per molecule. The finite-volume theory applies to two-point integrals over open subvolumes, and 
not to the KBIs computed from MDDFs.

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
    @test kbi(R; correction=:none, normalization=:bulk) ≈ R.kb
    @test kbi(R; correction=:G0, normalization=:bulk) == kbi(R; correction=:none, normalization=:bulk)
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
    g1 = kbi(R; correction=:G1, normalization=:bulk)
    g2 = kbi(R; correction=:G2, normalization=:bulk)
    w7 = kbi(R; correction=:W7, normalization=:bulk)
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

    # Not normalized results
    R = load("$data_dir/NAMD/water/rw_20_25.json")
    R.density.solvent_bulk = 0.0
    @test_throws "bulk solvent density is zero" kbi(R)
end
