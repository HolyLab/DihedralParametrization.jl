module DihedralParametrization

using BioStructures: Atom, Chain, Residue, atomname, bondangle, chain, collectatoms,
                     collectresidues, dihedralangle, omegaangles, phiangles, psiangles,
                     resid, residue, resname, sequentialresidues
using LinearAlgebra: cross, dot, norm
using StaticArrays: SMatrix, SVector

export atomcoordinates, bondparametrization, buildchain
export jacobianplan, coordinatejacobian, coordinatejacobian!, jtv!, jvp!, vhp!, vhp

include("tables.jl")
include("encode.jl")
include("decode.jl")
include("jacobian.jl")

end
