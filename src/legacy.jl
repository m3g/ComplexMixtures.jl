# some legacy interfaces, to not break old code

import Base: getproperty
getproperty(r::Result, s::Symbol) = getfield(r, Val(s))
getproperty(r::Result, ::Val{sum_md_count}) = r.coordination_number
getproperty(r::Result, ::Val{sum_md_count_random}) = r.coordination_number_random

const contrib = contributions

