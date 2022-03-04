"""

$(TYPEDEF)

Simple structure to contain the number of samples of each type
of calculation to compute final results

$(TYPEDFIELDS)

"""
@with_kw struct Samples
    md::Float64
    random::Int
end
