#
# ComplexMixtures.py
#
# A Python module to provide an inferface for the Julia ComplexMixtures.jl package.
# 
# See: https://m3g.github.com/ComplexMixtures.jl
#
# Author: L. Martinez / IQ-Unicamp, 2023.
#
import sys

#
# Check juliacall installation
#
try :
    print("Loading juliacall: on the first call this will install the julia executable.")
    from juliacall import Main as jl
except :
    print("""

    juliacall module not found. Install it with: 

    pip install juliacall

    """)
    sys.exit()

#
# Install and load necessary julia packages
#
try :
    jl.seval("import ComplexMixtures as cm")
    jl.seval("import PDBTools as pdb")
except :
    print("Installing the ComplexMixtures and PDBTools julia packages...")
    jl.Pkg.add("ComplexMixtures")
    jl.Pkg.add("PDBTools")
    jl.seval("import ComplexMixtures as cm")
    jl.seval("import PDBTools as pdb")

#
# Interfaces
#

# From PDBTools
readPDB = jl.pdb.readPDB
select = jl.pdb.select

# From ComplexMixtures
Selection = jl.cm.Selection
Trajectory = jl.cm.Trajectory
Options = jl.cm.Options
save = jl.cm.save
load = jl.cm.load
write = jl.cm.write
contrib = jl.cm.contrib
coordination_number = jl.cm.coordination_number
gr = jl.cm.gr
overview = jl.cm.overview

# For the possibly multi-threaded call to mddf, we need to disable the garbage collector
def mddf(*args, **kwargs) :
    jl.GC.enable(False)
    result = jl.cm.mddf(*args, **kwargs)
    jl.GC.enable(True)
    return result

# Covert python lists to julia arrays
def list(python_list) :
    return jl.map(jl.identity, python_list)















    
    

    
    


