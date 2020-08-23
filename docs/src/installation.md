# Installation

First you need to install the Julia language in your platform, from: 
[http://julialang.org](http://julialang.org). Julia version 1.3 or greater is required.

Next, run julia, and within the julia REPL interface, install the ComplexMixtures package using
```julia
julia> import Pkg

julia> Pkg.add("ComplexMixtures")

```
or simply
```julia
julia> ] add ComplexMixtures

```

The [PDBTools](http://m3g.iqm.unicamp.br/PDBTools) 
package will is one dependency that will be installed by default, and 
that will be used many times throughout the user guide. 

If you are first-time `julia` user, load these packages for the first
time after installation. Loading the `Plots` package, in particular, may
take quite a while when done for the first time, because it is compiled
at this point. To load the packages, use:

```
using ComplexMixtures, PDBTools, Plots
```

If no errors were shown in any of these steps, the packages are ready to
be used following the instructions and examples.

!!! tip
    The functions of the package are called, for example, using `ComplexMixtures.mddf(...)`.
    To avoid having to write `ComplexMixtures` all the time, define an
    accronym. For example:
    ```julia
    using ComplexMixtures ; const CM = ComplexMixtures
    CM.mddf(...)

    ```

