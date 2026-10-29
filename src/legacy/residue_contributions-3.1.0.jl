#=
    ResidueContributions_3_1_0

**LEGACY** support for reading `ResidueContributions` files saved with versions 2.8.0 - 3.1.0, 
in which the `version` field was named `Version`.

=#
@kwdef struct ResidueContributions_3_1_0
    Version::VersionNumber
    d::Vector{Float64}
    residue_contributions::Vector{Vector{Float64}}
    resnums::Vector{Int}
    xticks::Tuple{Vector{Int},Vector{String}}
end
_convert(rc::ResidueContributions_3_1_0) = ResidueContributions{MultipleResidueContribution}(
    d=rc.d,
    residue_contributions=rc.residue_contributions,
    resnums=rc.resnums,
    xticks=rc.xticks,
)
