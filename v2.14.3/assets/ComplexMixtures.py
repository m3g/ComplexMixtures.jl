#
# ComplexMixtures.py
#
# A Python module to provide an interface for the Julia ComplexMixtures.jl package.
# 
# See: https://m3g.github.com/ComplexMixtures.jl
#
# Author: L. Martinez / IQ-Unicamp, 2023.
#
# This script is adapted to version 2.13 of ComplexMixtures.jl
# and full functionality requires PDBTools.jl v3.0.0 or greater
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
read_pdb = jl.pdb.read_pdb
select = jl.pdb.select
select_with_vmd = jl.pdb.select_with_vmd
# Legacy compatibility
readPDB = jl.pdb.read_pdb

# From ComplexMixtures
AtomSelection = jl.cm.AtomSelection
Trajectory = jl.cm.Trajectory
SoluteGroup = jl.cm.SoluteGroup
SolventGroup = jl.cm.SolventGroup
Options = jl.cm.Options
save = jl.cm.save
load = jl.cm.load
write = jl.cm.write
contributions = jl.cm.contributions
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
