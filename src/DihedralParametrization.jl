module DihedralParametrization

using BioStructures: Atom, Chain, Residue, atomname, bondangle, chain, collectatoms,
                     collectresidues, dihedralangle, inscode, omegaangles, phiangles,
                     psiangles, resid, residue, resname, resnumber, sequentialresidues
using LinearAlgebra: cross, dot, norm
using StaticArrays: SMatrix, SVector

export BondParametrization, JacobianPlan
export atomcoordinates, atomcoordinates!, bondparametrization, buildchain, dihedralangles, dihedralangles!,
       dihedrallabels, natoms, ndihedrals
export jacobianplan, coordinatejacobian, coordinatejacobian!, vjp, vjp!, jvp, jvp!
export weightedhessian, weightedhessian!

VERSION >= v"1.11" && eval(Meta.parse("public AtomKey, DihedralLabel"))

include("tables.jl")
include("encode.jl")
include("decode.jl")
include("jacobian.jl")

end
