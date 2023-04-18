# Installation

## Install Julia

First you need to install the Julia language in your platform, from: 
[http://julialang.org](http://julialang.org). Julia version 1.6 or greater is required. Using the [juliaup](https://github.com/JuliaLang/juliaup) tool is also a highly recommended way of installing and keeping Julia up to date.

!!! tip
     New to Julia? Julia is a modern high-level yet performant programming language. Some tips
     and a nice workflow for using it effectively can be found [here](https://m3g.github.io/JuliaNotes.jl/stable/workflow/). 

     For this specific package, following a the step-by-step examples provided here after installing Julia should be enough. 

## Install the packages

Within Julia, to install the packages required for running the examples here you need to do:

```julia-repl
julia> import Pkg

julia> Pkg.add(["ComplexMixtures","Plots","PDBTools"])
```

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
julia> Pkg.add(["ComplexMixtures","PDBTools","Plots"])
```

That's it. Close Julia. Note that this created the files `Manifest.toml` and `Project.toml` in the `MyNewPaper` directory, which contain the information of packages and exact package versions you are using now on in this environment. Saving these files may be relevant for the future exact reproduction of your analyses. 

### Run your analysis scripts in that environment

Now, your analysis scripts, described in the next section in details, will look like: 

```julia
import Pkg; Pkg.activate("/home/user/Documents/MyNewPaper")

using ComplexMixtures
using PDBTools
using Plots

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
    To avoid having to write `ComplexMixtures` all the time, define an accronym. For example:
    ```julia
    import ComplexMixtures as CM
    CM.mddf(...)
    ```

