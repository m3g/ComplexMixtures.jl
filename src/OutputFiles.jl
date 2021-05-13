"""

$(TYPEDEF)

Structure to contain the names of the output files.

$(TYPEDFIELDS)

"""
@with_kw mutable struct OutputFiles

  output::String
  solute_atoms::String
  solvent_atoms::String

end

