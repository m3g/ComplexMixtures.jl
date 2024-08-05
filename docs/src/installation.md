# Installation

!!! note
    This is a package written in Julia. We invite you to experiment with the language, but if you want to just call this package 
    from Python, read the [From Python](@ref python) section of the manual. Understanding all the features of the package 
    requires reading the manual as whole. The syntaxes of using this package from Julia or Python are almost identical, and the
    motivation for using Python should be mostly the familiarity with further analysis tools, as the plotting packages. 

## Install Julia

First you need to install the Julia language, version 1.9 or greater is required. 
Using the [juliaup](https://github.com/JuliaLang/juliaup) tool is a highly recommended way of installing and keeping Julia up to date.

Alternatively, you can install Julia by downloading the binaries directly from [the Julia webpage](https://julialang.org).

!!! note
    New to Julia? Julia is a modern high-level yet performant programming language. Some tips
    and a nice workflow for using it effectively can be found [here](https://m3g.github.io/JuliaNotes.jl/stable/workflow/). 

    For this specific package, following a the step-by-step examples provided here after installing Julia should be enough. 

## Install the packages

Within Julia, to install the packages required for running the examples here you need to do:

```julia-repl
julia> import Pkg

julia> Pkg.add(["ComplexMixtures", "PDBTools", "Plots", "EasyFit", "LaTeXStrings"])
```

Here, [PDBTools.jl](https://m3g.github.io/PDBTools.jl) is an auxiliary package to read PDB files and select atoms within them.
The `Plots`, `EasyFit` and `LaTeXStrings` packages will help producing nice looking plots. 

Please read the recommended workflow below, for further information and to be sure to have a smoother experience.

## Recommended workflow for reproducibility

### Create an environment

Once Julia is installed, we recommend to create an environment that will contain all the packages you may use for your
analyses, including `ComplexMixtures`, in such a way that your results can always be reproduced and you don't get
any version incompatibility.

We illustrate this by creating the "MyNewPaper" environment, which will be hosted in a simple directory,
```bash
mkdir /home/user/Documents/MyNewPaper
```

Then, start Julia and activate the environment that will be hosted there:
```julia-repl
julia> import Pkg; Pkg.activate("/home/user/Documents/MyNewPaper")
  Activating new project at `~/Documents/MyNewPaper`
```

and add to this environment the packages that your analyses will require:

```julia-repl
julia> import Pkg; Pkg.add(["ComplexMixtures","PDBTools","Plots", "EasyFit", "LaTeXStrings"])
```

That's it. Close Julia. Note that this created the files `Manifest.toml` and `Project.toml` in the `MyNewPaper` directory, which contain the information of packages and exact package versions you are using now on in this environment. Saving these files may be relevant for the future exact reproduction of your analyses. 

### Run your analysis scripts in that environment

Now, your analysis scripts, described in the next section in details, will look like: 

```julia
import Pkg; Pkg.activate("/home/user/Documents/MyNewPaper")

using ComplexMixtures
using PDBTools
using Plots
using EasyFit
using LaTeXStrings

# etc ... 
```

And the script can be run with `julia -t auto script.jl` (where `-t auto` allows for multi-threading), or included in julia with `julia> include("./scritp.jl")`, as described in the next section.

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
    To avoid having to write `ComplexMixtures` all the time, define an acronym. For example:
    ```julia
    import ComplexMixtures as CM
    CM.mddf(...)
    ```

