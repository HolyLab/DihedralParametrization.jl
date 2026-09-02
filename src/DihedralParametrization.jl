"""
Represent protein chains by their rotatable-bond dihedral angles. Use
`bondparametrization` to construct the representation, `atomcoordinates` and
`dihedralangles` to convert between coordinates and angles, and
`jacobianplan` for derivatives with respect to the angles.

Numeric arrays representing coordinate lists, dihedrals, or derivatives must
use conventional one-based axes. Residue vectors and individual coordinate
vectors may use arbitrary axes.
"""
module DihedralParametrization

using BioStructures: AbstractAtom, AbstractResidue, Chain, atomname, atomnames, bondangle,
    collectatoms, collectresidues, coords, coords!, dihedralangle, inscode,
    omegaangles, phiangles, psiangles, resid, residue, resname, resnumber,
    sequentialresidues
using LinearAlgebra: cross, dot, norm
using StaticArrays: SMatrix, SVector

export BondParametrization, JacobianPlan, JacobianWorkspace
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
