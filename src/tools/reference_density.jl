"""
    ReferenceDensity

Structure that contains the estimates of the density of the solvent in different regions of 
the system, relative to the bulk density, as computed by [`reference_density`](@ref).

# Fields

- `d::Vector{Float64}`: distances to the solute (Å).
- `bulk::Float64`: the bulk density of the solvent, estimated beyond the cutoff (molecules/Å³).
- `outside::Vector{Float64}`: the density of the solvent beyond each distance, ``\\rho_{\\rm out}(d)``, 
   relative to `bulk`. This is the reference density of the Ganguly normalization of [`kbi`](@ref).
- `window::Vector{Float64}`: the density of the solvent between each distance and the cutoff,
   relative to `bulk`. It is not computed (`NaN`) for distances closer to the cutoff than the
   width of the shells, where the volume of the region is small and the estimate is noisy.
- `shell_d::Vector{Float64}` and `shell::Vector{Float64}`: the density of the solvent in shells
   of width `shell_width`, relative to `bulk`, and the distances of the centers of the shells.

!!! compat
    This structure was introduced in version 2.19.1.

"""
struct ReferenceDensity
    d::Vector{Float64}
    bulk::Float64
    outside::Vector{Float64}
    window::Vector{Float64}
    shell_d::Vector{Float64}
    shell::Vector{Float64}
end

function Base.show(io::IO, ::MIME"text/plain", rd::ReferenceDensity)
    i = findlast(isfinite, rd.window)
    print(io, chomp("""
    ReferenceDensity: densities relative to the bulk density ($(round(rd.bulk; sigdigits=6)) molecules/Å³)
        Density beyond d, at d = $(round(first(rd.d); digits=2)) Å: $(round(first(rd.outside); digits=5))
        Density between d = $(round(rd.d[i]; digits=2)) Å and the cutoff: $(round(rd.window[i]; digits=5))
    """))
end

"""
    reference_density(R::Result; shell_width::Real=1.0)

Compute the density of the solvent in different regions of the system, relative to the bulk density
of the solvent, `R.density.solvent_bulk`, which is estimated from the region beyond the cutoff. 
The result is a [`ReferenceDensity`](@ref) object.

For each distance ``d`` to the solute, the following densities are computed:

- `outside`: the density beyond ``d``, ``\\rho_{\\rm out}(d) = [N - N_{\\rm in}(d)]/[V - V(d)]``, 
  which is the reference density of the Ganguly normalization of [`kbi`](@ref).
- `window`: the density between ``d`` and the cutoff (not computed within `shell_width` of the cutoff). 
- `shell`: the density in shells of width `shell_width`.

The density beyond ``d`` is the reference density of the Ganguly normalization. Its deviation from the
bulk density is the deviation of the density between ``d`` and the cutoff, scaled by the fraction of the volume
beyond ``d`` that is within the cutoff. Thus, the density between ``d`` and the cutoff is the most sensitive 
indicator of the non-uniformity of the density, while the density beyond ``d`` is the one that affects the KBIs 
computed with the Ganguly normalization.

If the bulk density is properly estimated, all these ratios must be close to one, 
for ``d`` beyond the correlation length of the distribution. A density in the `window` that varies
as ``d`` approaches the cutoff, or that is systematically different from one, indicates that the
density of the solvent is not uniform beyond the correlation length, which may be caused by 
insufficient sampling or by long-range effects. In this case, the KBIs depend on the reference density
used for their normalization.

With `Plots` loaded, the result can be plotted with `plot(reference_density(R))`.

# Example

```julia-repl
julia> using ComplexMixtures, Plots

julia> R = load("result.json");

julia> rd = reference_density(R)

julia> plot(rd)
```

!!! compat
    This function was introduced in version 2.19.1.

"""
function reference_density(R::Result; shell_width::Real=1.0)
    _check_normalized(R)
    ρb = R.density.solvent_bulk
    Vd = cumsum(R.md_count_random) ./ ρb # Volume of the domain within d
    Nin = cumsum(R.md_count) # Number of solvent molecules within d
    outside = _ganguly_density(R) ./ ρb
    nb = max(1, round(Int, shell_width / R.files[1].options.binstep))
    window = fill(NaN, length(R.d))
    for i in 1:length(R.d)-nb
        ΔV = Vd[end] - Vd[i]
        ΔV > 0 && (window[i] = (Nin[end] - Nin[i]) / ΔV / ρb)
    end
    shell_d = Float64[]
    shell = Float64[]
    for k in 1:nb:length(R.d)-nb+1
        r = sum(@view(R.md_count_random[k:k+nb-1]))
        r > 0 || continue
        push!(shell_d, (R.d[k] + R.d[k+nb-1]) / 2)
        push!(shell, sum(@view(R.md_count[k:k+nb-1])) / r)
    end
    return ReferenceDensity(copy(R.d), ρb, outside, window, shell_d, shell)
end

@testitem "reference_density" begin
    using ComplexMixtures
    using ComplexMixtures: data_dir
    R = load("$data_dir/NAMD/water/rw_20_25.json")
    rd = reference_density(R)
    @test rd.bulk == R.density.solvent_bulk
    @test rd.outside ≈ ComplexMixtures._ganguly_density(R) ./ R.density.solvent_bulk
    i = findfirst(>=(20.0), R.d)
    ρb = R.density.solvent_bulk
    @test rd.window[i] ≈ sum(R.md_count[i+1:end]) / (sum(R.md_count_random[i+1:end]) / ρb) / ρb
    @test isnan(rd.window[end])
    @test all(isnan, rd.window[end-49:end])
    @test !isnan(rd.window[end-50])
    @test all(x -> isapprox(x, 1.0; atol=1e-3), filter(isfinite, rd.window[findfirst(>=(5.0), R.d):end-50]))
    @test all(x -> isapprox(x, 1.0; atol=0.05), rd.shell[findfirst(>=(5.0), rd.shell_d):end])
    @test length(rd.shell) == length(R.d) ÷ 50
    @test length(reference_density(R; shell_width=0.5).shell) == length(R.d) ÷ 25
    @test occursin("ReferenceDensity", sprint(show, MIME"text/plain"(), rd))
end
