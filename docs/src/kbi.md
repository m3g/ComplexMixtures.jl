```@meta
CollapsedDocStrings = true
```

# [Kirkwood-Buff integrals](@id kbi)

## Computing KBIs: the `kbi` function

The Kirkwood-Buff integrals (KBIs) of a `Result` object are obtained with the [`kbi`](@ref) function:

```julia
G = kbi(R)
```

which returns the KBI, in cm³ mol⁻¹, as a function of the upper limit of integration, ``L``, for 
each distance of `R.d`. The `R.kb` field of the `Result` contains the same values, `R.kb == kbi(R)`.
The `kbi` function allows choosing the corrections applied, and computes the KBIs with the 
current definitions for results saved with previous versions of the package, in which `R.kb` contained
the uncorrected KBI. By default, the KBI is computed with a weight function that corrects for the 
truncation of the integral at a finite distance, and with a reference density that corrects for 
the finite number of molecules in the simulation box, as explained below. The corrections apply 
to minimum-distance distribution functions (MDDFs), thus to solutes and solvents of any shape, and
to radial distribution functions (RDFs), which are MDDFs of single-atom solutes and solvents. 

The uncorrected (truncated) integral of the MDDF is obtained with `kbi(R; correction=:none, normalization=:mddf)`.

!!! compat
    The `kbi` function was introduced in version 2.19.0. Its application to MDDFs, the `:W7` correction,
    and the `normalization` option were introduced in version 2.19.1, in which the defaults were set 
    to `correction=:W7` and `normalization=:ganguly`, and `R.kb` was set to `kbi(R)`. In previous versions,
    `R.kb` contained the uncorrected KBI.

## Theory in brief

### The truncated KBI

The KBI of a pair of species is, in the thermodynamic limit, the integral of the excess density of 
the solvent around the solute over all space. In terms of the minimum-distance distribution 
function, ``g(d)``, 

```math
G_\infty = \int_0^\infty \left[g(d) - 1\right] \frac{dV(d)}{dd}\, dd
```

where ``V(d)`` is the volume of the domain within a minimum distance ``d`` of the solute (for 
single-atom solutes, ``dV/dd = 4\pi d^2``, and ``g(d)`` is the RDF). In a simulation, ``g(d)`` is 
known only up to a finite distance, and the integral must be truncated at ``L``. The truncated 
integral, ``G_0(L)``, converges poorly with ``L``: the oscillations of ``g(d) - 1`` are amplified 
by the volume element, and the integral oscillates around its limiting value even at long distances.

### Weight functions

Krüger and Vlugt [1,2] derived, from the theory of finite-volume KBIs, weight functions that estimate
``G_\infty`` from the RDF known up to ``L``:

```math
G_W(L) = \int_0^L \left[g(d) - 1\right] \frac{dV(d)}{dd}\, W(d/L)\, dd
```

Santos [3] showed that these weights follow from a purely mathematical identity, valid for any 
one-dimensional integral, and generalized them to a family of weights, ``W_n^{(k)}(x)``, that vanish
more smoothly at ``x = d/L = 1``. Since the identity does not depend on the origin of the integrand, 
the weights can be applied to the KBI computed from MDDFs. The available weights are:

| `correction` | Weight, ``W(x)``, ``x = d/L`` | |
|:---------:|:-------------------------------|:---------|
| `:none` (or `:G0`) | ``1`` | The truncated integral. Large oscillations. |
| `:G1` | ``1 - x^3`` | Ref. [1]. |
| `:G2` | ``1 - \frac{23}{8}x^3 + \frac{3}{4}x^4 + \frac{9}{8}x^5`` | Ref. [2], Eq. 24. |
| `:W7` | ``(1-x)^4 (1 + \frac{35}{16}x)(1 + \frac{1225}{256}x^2)(1 + \frac{29}{16}x + \frac{5}{4}x^2 + \frac{5}{16}x^3)`` | Ref. [3]. Default. |

The weights go to zero at ``d = L``, thus the oscillations of ``g(d)`` near the truncation point do
not propagate to the integral:

```@example kbi
using Plots
x = range(0, 1, length=200)
W7(x) = (1-x)^4 * (1 + 35x/16) * (1 + 1225x^2/256) * (1 + 29x/16 + 5x^2/4 + 5x^3/16)
plot(x, ones(length(x)); label="G₀ (:none)", linewidth=2)
plot!(x, @. 1 - x^3; label="G₁ (:G1)", linewidth=2)
plot!(x, @. 1 - 23/8 * x^3 + 3/4 * x^4 + 9/8 * x^5; label="G₂ (:G2)", linewidth=2)
plot!(x, W7.(x); label="W₇⁽³⁾ (:W7)", linewidth=2)
plot!(xlabel="x = d/L", ylabel="W(x)", framestyle=:box, size=(500, 350))
```

The weights correct the truncation of the integral when ``g(d) - 1`` oscillates around zero at ``L``. 
They do not correct slowly decaying (monotonic) tails of the distribution, nor errors in the reference density.

### Reference density

The distribution function is the ratio between the density of the solvent at each distance and a 
reference density, which must be the density of the solvent in the bulk. In a simulation with a fixed number
of molecules, the accumulation (or depletion) of solvent molecules around the solute changes the
density of the rest of the box, and an inaccurate reference density causes a drift of the KBI at long 
distances, because the volume element grows with ``L``. Two normalizations are available: 

- `normalization=:mddf`: the reference density is the bulk density of the solvent, `R.density.solvent_bulk`, 
  estimated from the region beyond the cutoff, at all distances. This is the normalization of the MDDF,
  `R.mddf`, thus the KBI is the integral of `R.mddf`.
- `normalization=:ganguly` (default): the reference density at each distance ``d`` is the density 
  of the solvent outside the domain within ``d``,
  ```math
  \rho_{\rm ref}(d) = \frac{N - N_{\rm in}(d)}{V - V(d)}
  ```
  where ``N`` is the number of solvent molecules (minus one if the solute and the solvent are the same), 
  ``N_{\rm in}(d)`` is the average number of solvent molecules within ``d``, and ``V`` is the volume of the box.
  This is the correction for closed systems proposed by Ganguly and van der Vegt [4].

If the two normalizations give different KBIs at the distances of interest, the estimate of the 
reference density is a source of error that must be considered. The density of the solvent in different 
regions of the system can be inspected with the [`reference_density`](@ref) function (see 
[Inspecting the reference density](@ref kbi_reference_density) below).

## Example: water

Here we use the distribution functions of a simulation of pure water (6845 molecules, cubic box of 
~58.7 Å), computed up to 25 Å. In `rw_20_25.json` the solute and the solvent are the complete water 
molecules, thus the distribution is an MDDF. In `rwO_20_25.json` they are the oxygen atoms only,
thus the distribution is the O–O RDF:

```@example kbi
using ComplexMixtures
using ComplexMixtures: data_dir
R = load(joinpath(data_dir, "NAMD", "water", "rw_20_25.json"))
RO = load(joinpath(data_dir, "NAMD", "water", "rwO_20_25.json"))
R.solute.natomspermol, RO.solute.natomspermol
```

### Running KBIs

The [`kbi`](@ref) function returns the KBIs as a function of the upper limit of integration, `L`, 
which corresponds to the distances of `R.d`:

```@example kbi
g0 = kbi(R; correction=:none, normalization=:mddf) # truncated
g2 = kbi(R; correction=:G2)
w7 = kbi(R) # correction=:W7, normalization=:ganguly
w7O = kbi(RO) # O-O RDF
plot(R.d, g0; label="G₀ (truncated)", linewidth=2, color=:gray)
plot!(R.d, g2; label="G₂", linewidth=2)
plot!(R.d, w7; label="W₇⁽³⁾ (default)", linewidth=2, color=:black)
plot!(RO.d, w7O; label="W₇⁽³⁾, O–O RDF", linewidth=2, linestyle=:dash)
plot!(xlabel="L / Å", ylabel="KBI / cm³ mol⁻¹", ylims=(-30, 0), xlims=(3, 25), framestyle=:box, size=(600, 400))
```

The truncated integral, ``G_0``, oscillates significantly, and drifts upwards at long distances. The 
weighted integrals reach a plateau at about 8 Å. A good estimate of ``G_\infty`` is the average value 
of the weighted KBI in the plateau:

```@example kbi
using Statistics
plateau = findall(L -> 8 <= L <= 15, R.d)
mean(w7[plateau]), std(w7[plateau])
```

The value obtained from the MDDF is similar to that obtained from the O–O RDF:

```@example kbi
plateau = findall(L -> 10 <= L <= 20, RO.d)
mean(w7O[plateau]), std(w7O[plateau])
```

At distances larger than ~15 Å, the KBI computed from the MDDF drifts upwards. The weights and the Ganguly 
normalization reduce the drift, but do not remove it. Since this drift is not a truncation effect, it is 
probably associated with the sampling or with the estimate of the reference density, and the values at 
shorter distances are preferable. 

As a reference, the KBI of a pure liquid is related to its isothermal compressibility, ``\kappa_T``, by 
``G_\infty = RT\kappa_T - 1/\rho``, where ``\rho`` is the molar density. For water at 298 K, the experimental 
value is ``G_\infty \approx 1.1 - 18.1 \approx -17`` cm³ mol⁻¹. The exact value obtained from the simulation 
depends on the water model, and on the sampling.

### Finite-volume KBIs and extrapolation

For radial distribution functions, an alternative estimate of ``G_\infty`` is provided by the theory
of finite-volume KBIs [1,2]. The KBI of a sphere of diameter ``L`` is

```math
G(L) = \int_0^L h(r)\, 4\pi r^2 \left(1 - \frac{3}{2}x + \frac{1}{2}x^3\right) dr, \quad x = \frac{r}{L},
```

where ``r`` is the distance between two points *inside* the sphere. For large ``L``, 
``G(L) = G_\infty + F_\infty/L + O(1/L^2)``, where ``F_\infty`` is a surface term. ``G(L)`` is not an estimate 
of ``G_\infty``, but ``G_\infty`` can be obtained by extrapolating ``G(L)`` to ``1/L \to 0``. 

This theory applies to two-point integrals over open subvolumes, and not to the KBIs computed from MDDFs,
which are integrals around a fixed solute. Thus, [`finite_volume_kbi`](@ref) requires a single-atom solute
(it issues a warning otherwise), and uses the distribution of the distances to the reference atom of the solvent
(`R.rdf_count`).

The [`finite_volume_kbi`](@ref) function computes the finite-volume KBIs of spheres of diameter ``L``,
for all ``L`` up to the cutoff. The plot of ``G(L)`` as a function of ``1/L`` must be linear for large enough ``L``:

```@example kbi
fv = finite_volume_kbi(RO)
```

```@example kbi
plot(1 ./ fv.L, fv.G; label="G(L)", linewidth=2)
plot!(xlabel="1/L / Å⁻¹", ylabel="G(L) / cm³ mol⁻¹", xlims=(0, 0.25), ylims=(-20, 0), framestyle=:box, size=(600, 400))
```

The function [`extrapolate_kbi`](@ref) fits the linear relation to the data in a range of values of ``L``,
and returns the intercept, ``G_\infty``, and the slope, ``F_\infty``. Here, we use ``10 \leq L \leq 20`` Å:

```@example kbi
ext = extrapolate_kbi(fv, (10.0, 20.0))
```

The line obtained can be compared with the data, to check if the relation is linear in the range chosen:

```@example kbi
x = range(0, 0.1, length=10)
plot!(x, ext.Ginf .+ ext.F .* x; label="Fit, 10 ≤ L ≤ 20 Å", linewidth=2, linestyle=:dash)
vspan!([1 / 20, 1 / 10]; alpha=0.2, color=:gray, label=nothing)
```

The extrapolated value is consistent with the plateau of the weighted KBIs. The result of the extrapolation should 
not depend much on the range of ``L`` used. This must be checked, because the linear relation is valid
only for ``L`` larger than the correlation length of the fluid, and the data at large ``L`` may be 
affected by the drift discussed above:

```@example kbi
for Lrange in ((8.0, 16.0), (10.0, 20.0), (10.0, 25.0), (15.0, 25.0))
    (; Ginf, F) = extrapolate_kbi(fv, Lrange)
    println("L range: $Lrange Å: Ginf = $(round(Ginf; digits=2)) cm³/mol, F = $(round(F; digits=1)) cm³ Å/mol")
end
```

## [Inspecting the reference density](@id kbi_reference_density)

The [`reference_density`](@ref) function computes the density of the solvent in different regions of the
system, relative to the bulk density (estimated beyond the cutoff): the density in shells around the solute, 
the density between each distance ``d`` and the cutoff, and the density beyond ``d``, which is the reference
density of the Ganguly normalization:

```@example kbi
rd = reference_density(R)
```

```@example kbi
plot(rd; size=(600, 400))
```

The plot shows the deviations of these densities from the bulk density, in percent. Beyond the correlation 
length of the distribution, all deviations should be close to zero. Here, for water, they are smaller than 0.1%,
except for the noise of the density in the shells. The density beyond ``d`` is the reference density of the
Ganguly normalization, and its deviation is the one that affects the KBIs. The density between ``d`` and the
cutoff is a more sensitive indicator of the non-uniformity of the density, because the deviation of the density
beyond ``d`` is diluted by the volume beyond the cutoff. If the
density between ``d`` and the cutoff varies as ``d`` approaches the cutoff, or if the densities are systematically
different from one, the density of the solvent is not uniform beyond the correlation length. This may be caused 
by insufficient sampling or by long-range effects, and the KBIs will depend on the reference density used to
normalize them.

## Interpreting the results

- **Use `kbi(R)`.** Plot the KBI as a function of ``L``, and report the value (or average) in the 
  range where it is stable. Compare the results obtained with the different corrections and normalizations:
  the differences between them are an indication of the uncertainty of the result.

- **Drifts at long distances.** If all estimates drift at long distances, the drift is not a truncation
  effect. It usually indicates an imprecise estimate of the reference density, insufficient sampling, or 
  a slowly decaying tail of the distribution. In this case, use the values at shorter distances, where
  the estimates are stable, increase the sampling (for example, by performing multiple independent 
  simulations and merging the results with [`merge`](@ref)), or increase the size of the simulation box.

- **Extrapolation.** For RDFs, the extrapolation of the finite-volume KBIs is a cross-check of the 
  weighted estimates. It relies on a fit whose result depends on the range of ``L``, and involves an 
  extrapolation from ``1/L \sim 0.05`` Å⁻¹ to zero, which amplifies the noise.

- **Group contributions.** The decomposition of the KBI into contributions of groups of atoms
  ([`contributions`](@ref) with `type=:kbi`) accepts the same `correction` and `normalization` options,
  and the contributions sum to the KBI computed by `kbi` with the same options.

## Reference functions

```@docs
kbi
reference_density
ComplexMixtures.ReferenceDensity
Plots.plot(::ComplexMixtures.ReferenceDensity)
finite_volume_kbi
extrapolate_kbi
ComplexMixtures.FiniteVolumeKBI
```

## References

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
