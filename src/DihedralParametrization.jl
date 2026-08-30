module DihedralParametrization

using BioStructures
using OrderedCollections
using LinearAlgebra
using StaticArrays

export atomcoordinates, bondparametrization, buildchain
export jacobianplan, coordinatejacobian, coordinatejacobian!, jtv!, jvp!

include("tables.jl")
include("encode.jl")
include("decode.jl")
include("jacobian.jl")

end
