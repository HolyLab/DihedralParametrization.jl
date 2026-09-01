```@meta
CurrentModule = DihedralParametrization
```

# Derivatives

For `X = atomcoordinates(bp, θ, n, cα, c)`, the package provides the
Jacobian `J = ∂X/∂θ` and products involving its first and second
derivatives. Each rotatable dihedral rotates its dependent atoms about a fixed
axis. Thus `∂X[t]/∂θ_k` is zero when `k` does not affect atom `t`; otherwise,
it is the cross product of the axis direction and the displacement from the
axis to `X[t]`.

## `jacobianplan`

[`jacobianplan`](@ref) records the dependency tree and rotation axes from
`bp`'s build sequence. The resulting [`JacobianPlan`](@ref) can be reused for
any configuration built from the same `bp`.

`jacobianplan` rejects build steps whose dependencies cannot be represented
as nested rigid rotations.

## Jacobian, and Jacobian/vector products

- [`coordinatejacobian`](@ref) / [`coordinatejacobian!`](@ref) build the
  dense `3*natoms × ndih` matrix `J = ∂X/∂θ`, evaluated at a specific `X`.
- [`jvp!`](@ref) computes `δx = J * v` (a Jacobian-vector product) without
  forming `J`.
- [`jtv!`](@ref) computes `g = J' * w` (a vector-Jacobian product) without
  forming `J`.
- [`vhp!`](@ref) / [`vhp`](@ref) compute
  `S[i,j] = Σ_a w[a] ⋅ ∂²X[a]/∂θ_i∂θ_j`, the contraction of the coordinate
  map's second derivative with a per-atom cotangent `w`, without forming the
  second-derivative tensor.

These functions take a `JacobianPlan` and coordinates `X`; call
`atomcoordinates` first. The vector products also take a per-atom vector `w`
or per-dihedral vector `v`.

## Worked example: gradient of a coordinate objective

This uses `jtv!` to differentiate a squared-distance objective:

```julia
using DihedralParametrization, BioStructures

struc = read("AF-M3YHX5-F1-model_v4_hydrogens.cif", MMCIFFormat)
specializeresnames!(struc)
chain = struc["A"]
bp, dihedrals = bondparametrization(chain)
plan = jacobianplan(bp)

n, cα, c = chain[1]["N"].coords, chain[1]["CA"].coords, chain[1]["C"].coords

function objective_and_grad(bp, plan, θ, n, cα, c, target)
    X = atomcoordinates(bp, θ, n, cα, c)
    w = X .- target                     # gradient of (1/2) Σ |X - target|² w.r.t. X
    val = 0.5 * sum(v -> sum(abs2, v), w)
    g = jtv!(zeros(plan.ndih), plan, X, w)
    return val, g
end

target = atomcoordinates(bp, dihedrals, n, cα, c)  # a structure to match
val, g = objective_and_grad(bp, plan, dihedrals, n, cα, c, target)
```

`g` is the gradient of `val` with respect to `dihedrals`.

## AD interop

`atomcoordinates` also works with automatic differentiation through
[DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl)
and a backend such as [Mooncake.jl](https://github.com/compintell/Mooncake.jl)
or [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl). These
backends differentiate the build steps directly and do not use
`jacobianplan`:

```julia
using DihedralParametrization, DifferentiationInterface
import Mooncake

bp, dihedrals = bondparametrization(chain)
flatten(X) = reduce(vcat, Vector.(X))
f(θ) = flatten(atomcoordinates(bp, θ, chain))

backend = DifferentiationInterface.AutoMooncakeForward(; config=nothing)
J = jacobian(f, backend, dihedrals)   # 3*natoms × ndih, comparable to coordinatejacobian(plan, X)
```

`ForwardDiff.jacobian(f, dihedrals)` works with the same `f`.

## Conventions

- **Units.** All angles are in radians: the dihedrals passed to
  `atomcoordinates`, the fixed values stored in `bp`, and the `θbcd`/`ϕbc`
  arguments of `snnerf`.

- **Layout and ordering of `dihedrals`.** For a chain of `nres` residues,
  `dihedrals` (of length `ndihedrals(bp)`) is built in two passes:

  1. **Backbone**, residue by residue from residue 1 to residue `nres-1`: for
     each residue `i`, the entry for ψ_i (the dihedral about the Cα_i–C_i
     bond, which places N of residue `i+1`) always appears, immediately
     followed by the entry for φ_{i+1} (about the N_{i+1}–Cα_{i+1} bond,
     which places C of residue `i+1`) *only if* `bp.phirotatable[i+1]` is
     `true`. After all `nres-1` residues, one final entry places the
     terminal OXT atom (about the Cα_nres–C_nres bond) and always appears.
  2. **Sidechains**, residue by residue from residue 1 to residue `nres`: for
     each residue, one entry per rotatable `Extend` build step in that
     residue's build sequence, in the sequence's own order (typically χ1,
     χ2, … outward from Cα, though a residue's sequence can also expose a
     single-bond rotation that is not conventionally numbered as a χ angle,
     such as a terminal methyl's orientation).

- **Fixed vs. rotatable dihedrals.** Three kinds of dihedral angle are fixed
  by `bp` and never appear in `dihedrals`:

  + Every ω (the dihedral about each peptide C–N bond), stored in
    `bp.omegas`.
  + φ for a residue whose ring fixes it, such as proline: there
    `bp.phirotatable[i]` is `false` and the fixed value is `bp.phi[i]`.
    Proline's ring also fixes its χ1 and χ2, which appear as non-rotatable
    `Extend` steps rather than entries in `dihedrals`.
  + Any other `Extend` build step marked non-rotatable, whose fixed value is
    stored in that step's own `ϕ` field.
