import Pkg
Pkg.activate(temp=true)
Pkg.develop("ComplexMixtures")
Pkg.add("WGLMakie")
Pkg.add("Bonito")
using ComplexMixtures, WGLMakie, Bonito
dir="/home/leandro/Documents/ComplexMixturesExamples/scripts/example1"
ComplexMixtures.gui(;
    pdbfile="$dir/system.pdb",
    result="$dir/glyc50_results.json"
)
