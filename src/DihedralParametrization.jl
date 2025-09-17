module DihedralParametrization

using BioStructures
using OrderedCollections
using LinearAlgebra
using StaticArrays

export atomcoordinates, bondparametrization, buildchain

include("tables.jl")
include("encode.jl")
include("decode.jl")

end
