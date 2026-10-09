```@meta
CollapsedDocStrings = true
```

# [Kirkwood-Buff integrals: convergence and finite-size corrections](@id kbi)

The Kirkwood-Buff integral (KBI) of a pair of species is, in the thermodynamic limit, the integral of the 
excess density of the solvent around the solute, over all space:

```math
G_\infty = \int_0^\infty h(r)\, 4\pi r^2\, dr, \quad h(r) = g(r) - 1.
```

In a simulation, ``g(r)`` is known only up to a finite distance ``L`` (the `cutoff` of the calculation), 
and the integral must be truncated. The `Result` object returned by `mddf` contains two such 
truncated integrals, as a function of the distance:

- `R.kb`: the KBI computed from the minimum-distance distribution function (MDDF).
- `R.kb_rdf`: the KBI computed from the distribution of the distances to the reference atom of the solvent. 
  If the solute has a single atom per molecule, this is the radial distribution function (RDF) of the 
  pair, and `R.kb_rdf` is the usual "running" KBI.

The running KBI converges poorly with ``L``: the oscillations of ``h(r)`` are amplified by the 
``r^2`` factor, and the integral oscillates around its limiting value even at long distances. 
For radial distribution functions, the convergence can be greatly improved by using the 
results of the theory of finite-volume KBIs, as proposed by Krüger and Vlugt [1,2], 
implemented in the [`kbi`](@ref), [`finite_volume_kbi`](@ref) and [`extrapolate_kbi`](@ref) functions,
which are described here.

!!! compat
    The functions described in this page were introduced in version 2.19.0.

## Theory in brief

### Finite-volume KBIs

The KBI of a finite (sub)volume ``V``, ``G(V)``, measures the particle number fluctuations
inside ``V``. It can be written as a radial integral with a purely geometrical weight, which, for a
sphere of diameter ``L``, is [1,2]:

```math
G(L) = \int_0^L h(r)\, 4\pi r^2 \left(1 - \frac{3}{2}x + \frac{1}{2}x^3\right) dr, \quad x = \frac{r}{L}.
```

Note that ``r`` is the distance between two points *inside* the sphere, thus it varies from ``0`` to the 
*diameter* ``L``, and ``G(L)`` requires ``g(r)`` up to ``r = L``. For large ``L``, [2]

```math
G(L) = G_\infty + \frac{F_\infty}{L} + O\left(\frac{1}{L^2}\right),
```

where ``F_\infty`` is a surface term. ``G(L)`` is a smooth function of ``L``, but it is *not*
an estimate of ``G_\infty``: it differs from it by ``F_\infty/L``, which decays slowly. ``G_\infty`` can, 
however, be obtained by extrapolating ``G(L)`` as a function of ``1/L`` to ``1/L \to 0``.

### Improved estimators of the infinite-volume KBI

Using the relation above, Krüger and Vlugt derived weight functions that estimate ``G_\infty`` directly
from ``g(r)`` known up to ``L`` [1,2]: 

```math
G_k(L) = \int_0^L h(r)\, u_k(r)\, dr
```

with

| Estimator | Weight, ``u_k(r)``, ``x = r/L`` | Behavior |
|:---------:|:-------------------------------|:---------|
| ``G_0`` | ``4\pi r^2`` | The truncated integral. Large oscillations. |
| ``G_1`` | ``4\pi r^2 (1 - x^3)`` | Ref. [1]. Smaller oscillations. |
| ``G_2`` | ``4\pi r^2 \left(1 - \frac{23}{8}x^3 + \frac{3}{4}x^4 + \frac{9}{8}x^5\right)`` | Ref. [2], Eq. 24. Error ``\sim 1/L^3``. Recommended. |

The weights of ``G_1`` and ``G_2`` go to zero at ``r = L``, smoothly in the case of ``G_2``. Thus, the 
oscillations of ``h(r)`` near the truncation point do not propagate to the integral, which is the 
reason of the poor convergence of ``G_0``. The weight functions are shown below:

```@example kbi
using Plots
x = range(0, 1, length=200)
plot(x, ones(length(x)); label="G₀", linewidth=2)
plot!(x, @. 1 - x^3; label="G₁", linewidth=2)
plot!(x, @. 1 - 23/8 * x^3 + 3/4 * x^4 + 9/8 * x^5; label="G₂", linewidth=2)
plot!(x, @. 1 - 3/2 * x + 1/2 * x^3; label="Finite-volume sphere, G(L)", linewidth=2, linestyle=:dash)
plot!(xlabel="x = r/L", ylabel="u(r) / 4πr²", framestyle=:box, size=(500, 350))
```

### When can these corrections be used

- **Radial distribution functions only.** The theory applies to functions of the distance between two
  points (two atoms). Thus, the solute must be defined with a single atom per molecule (for example, the 
  oxygen atom of water). The solvent may contain more than one atom: the distances are computed to its 
  reference atom (by default, the first atom of the molecule). If the solute has more than one atom per 
  molecule, the functions will issue a warning, because the distribution is a minimum-distance
  distribution, for which the theory does not apply.

- **Truncation, not ensemble, errors.** The corrections address the truncation of the integral at a 
  finite distance. They do not correct the systematic error of ``g(r)`` computed from simulations of
  closed systems (with a fixed number of molecules), nor errors in the estimate of the bulk density. 
  These errors are small at each distance, but are amplified by the ``r^2`` factor, and appear as a 
  drift of all estimators at long distances. The values should be taken from a range of distances 
  where the estimates are stable.

## Example: water

Here we use the oxygen-oxygen distribution function of a simulation of pure water (6845 molecules,
cubic box of ~58.7 Å), computed up to 25 Å. The solute and the solvent were defined
by the water oxygen atoms only, thus the distribution is a radial distribution function:

```@example kbi
using ComplexMixtures
using ComplexMixtures: data_dir
R = load(joinpath(data_dir, "NAMD", "water", "rwO_20_25.json"))
R.solute.natomspermol # must be one
```

### Improved running KBIs

The [`kbi`](@ref) function returns the running KBIs computed with the different estimators, as a function
of the upper limit of integration, `L`, which corresponds to the distances of `R.d`:

```@example kbi
g0 = kbi(R; correction=:G0) # same as R.kb_rdf
g1 = kbi(R; correction=:G1)
g2 = kbi(R; correction=:G2)
plot(R.d, g0; label="G₀ (truncated)", linewidth=2)
plot!(R.d, g1; label="G₁", linewidth=2)
plot!(R.d, g2; label="G₂", linewidth=2)
plot!(xlabel="L / Å", ylabel="KBI / cm³ mol⁻¹", ylims=(-30, 0), xlims=(3, 25), framestyle=:box, size=(600, 400))
```

The truncated integral, ``G_0``, oscillates significantly, and its value depends strongly on the distance 
at which it is read. ``G_1`` and, particularly, ``G_2`` are much more stable, and reach a plateau 
already at about 8 Å. A good estimate of ``G_\infty`` is the average value of ``G_2`` in the plateau:

```@example kbi
using Statistics
plateau = findall(L -> 10 <= L <= 20, R.d)
mean(g2[plateau]), std(g2[plateau])
```

At distances larger than ~18 Å, an upward drift is visible in ``G_0`` and, to a lesser extent, in ``G_1``, 
while ``G_2`` is barely affected. Since this drift is not a truncation effect, it is probably associated
with the closed-system errors mentioned above, and the estimates at shorter distances are preferable.

As a reference, the KBI of a pure liquid is related to its isothermal compressibility, ``\kappa_T``, by 
``G_\infty = RT\kappa_T - 1/\rho``, where ``\rho`` is the molar density. For water at 298 K, the experimental 
value is ``G_\infty \approx 1.1 - 18.1 \approx -17`` cm³ mol⁻¹. The exact value obtained from the simulation 
depends on the water model, and on the sampling.

### Finite-volume KBIs and extrapolation

The [`finite_volume_kbi`](@ref) function computes the finite-volume KBIs of spheres of diameter ``L``,
for all ``L`` up to the cutoff. Since ``G(L) = G_\infty + F_\infty/L``, the plot of ``G(L)`` as a function 
of ``1/L`` must be linear for large enough ``L``:

```@example kbi
fv = finite_volume_kbi(R)
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

The extrapolated value is consistent with the plateau of ``G_2``. The result of the extrapolation should 
not depend much on the range of ``L`` used. This must be checked, because the linear relation is valid
only for ``L`` larger than the correlation length of the fluid, and the data at large ``L`` may be 
affected by the drift discussed above:

```@example kbi
for Lrange in ((8.0, 16.0), (10.0, 20.0), (10.0, 25.0), (15.0, 25.0))
    (; Ginf, F) = extrapolate_kbi(fv, Lrange)
    println("L range: $Lrange Å: Ginf = $(round(Ginf; digits=2)) cm³/mol, F = $(round(F; digits=1)) cm³ Å/mol")
end
```

## Interpreting the results

- **Prefer ``G_2``.** The `:G2` estimator of `kbi` is the simplest and most robust estimate of ``G_\infty``. 
  Plot it as a function of ``L``, and report the value (or average) in the range where it is stable.

- **Use the extrapolation as a cross-check.** The extrapolation of ``G(L)`` relies on a fit, whose result 
  depends on the range of ``L``, and involves an extrapolation from ``1/L \sim 0.05`` Å⁻¹ to zero, 
  which amplifies the noise. When ``G_2`` and the extrapolation agree, the result is reliable. 

- **The surface term, ``F_\infty``**, is a property of the fluid, and describes how the density fluctuations 
  in a finite volume deviate from those of the infinite system. It has units of cm³ mol⁻¹ Å here. 

- **Drifts at long distances.** If all estimators drift at long distances, the drift is not a truncation
  effect, and these corrections cannot remove it. It usually indicates a closed-system error in ``g(r)``,
  or an imprecise estimate of the bulk density. In this case, use the values at shorter distances, where
  the estimates are stable, increase the size of the simulation box, or consider corrections for closed
  systems [3].

- **Minimum-distance distributions.** For solutes with more than one atom (proteins, polymers, etc.), use the KBI
  computed from the MDDF, `R.kb`. The corrections described here do not apply. The convergence of MDDF-based 
  KBIs is discussed in the [Kirkwood-Buff integrals and convergence](@ref) section.

## Reference functions

```@docs
kbi
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
3. P. Ganguly, N. F. A. van der Vegt, Convergence of Sampling Kirkwood–Buff Integrals of Aqueous Solutions 
   with Molecular Dynamics Simulations. *J. Chem. Theory Comput.* 9, 1347 (2013).
   [DOI: 10.1021/ct301017q](https://doi.org/10.1021/ct301017q)
