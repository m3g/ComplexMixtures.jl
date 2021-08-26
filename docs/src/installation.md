# Installation

First you need to install the Julia language in your platform, from: 
[http://julialang.org](http://julialang.org). Julia version 1.6 or greater is required.

Next, run julia, and within the julia REPL interface, install the ComplexMixtures package using
```julia
julia> import Pkg

julia> Pkg.add("ComplexMixtures")

```
or simply
```julia
julia> ] add ComplexMixtures

```

To follow all the examples provided in this manual, the 
[PDBTools](http://m3g.iqm.unicamp.br/PDBTools) 
and [Plots](http://docs.juliaplots.org/latest/) have to be installed as well:
```julia
julia> ] add PDBTools, Plots

```

If you are first-time `julia` user, load these packages for the first
time after installation. Loading the `Plots` package, in particular, may
take quite a while when done for the first time, because it is compiled
at this point (this was greatly improved in Julia versions greater than 1.6, which
are highly recommended). To load the packages, use:

```
using ComplexMixtures, PDBTools, Plots
```

If no errors were shown in any of these steps, the packages are ready to
be used following the instructions and examples.

!!! tip
    By loading the package with 
    ```julia
    using ComplexMixtures

    ```
    the most common functions of the package become readily available by their direct name, 
    for example `mddf(...)`.

    If you don't want to bring the functions into the scope of your script, use
    ```julia
    import ComplexMixtures

    ```
    Then, the functions of the package are called, for example, using `ComplexMixtures.mddf(...)`.
    To avoid having to write `ComplexMixtures` all the time, define an
    accronym. For example:
    ```julia
    import ComplexMixtures 
    const CM = ComplexMixtures
    CM.mddf(...)

    ```

