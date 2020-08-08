# Installation

First you need to install the Julia language in your platform, from: 
[http://julialang.org](http://julialang.org)

Next, run julia, and within the julia REPL interface, install the MDDF package using
```julia
julia> import Pkg

julia> Pkg.add("https://github.com/m3g/MDDF")

```
The [PDBTools](http://m3g.iqm.unicamp.br/PDBTools) 
package will is one dependency that will be installed by default, and 
that will be used many times throughout the user guide. 

If you are first-time `julia` user, load these packages for the first
time after installation. Loading the `Plots` package, in particular, may
take quite a while when done for the first time, because it is compiled
at this point. To load the packages, use:

```
using MDDF, PDBTools, Plots
```

If no errors were shown in any of these steps, the packages are ready to
be used following the instructions and examples.
