```@meta
CurrentModule = DihedralParametrization
```

# Derivatives

For `X = atomcoordinates(bp, θ, (n, cα, c))`, the package provides the
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
- [`jvp`](@ref) / [`jvp!`](@ref) compute `δx = J * v` (a Jacobian-vector
  product) without forming `J`.
- [`vjp`](@ref) / [`vjp!`](@ref) compute `g = J' * w` (a vector-Jacobian
  product) without forming `J`.
- [`weightedhessian`](@ref) / [`weightedhessian!`](@ref) compute
  `S[i,j] = Σ_a w[a] ⋅ ∂²X[a]/∂θ_i∂θ_j`, the Hessian with respect to θ of
  the scalar `Σ_a w[a] ⋅ X[a](θ)`, evaluated at `X`. It is the full
  `ndih × ndih` Hessian of `w ⋅ X`, not a Hessian-vector product.

These functions take a `JacobianPlan` and coordinates `X`; call
`atomcoordinates` (or, to reuse a buffer, `atomcoordinates!`) first. The
vector products also take a per-atom vector `w` or per-dihedral vector `v`.

## Worked example: gradient of a coordinate objective

This uses `vjp!` to differentiate a squared-distance objective:

```julia
using DihedralParametrization, BioStructures

struc = read("AF-M3YHX5-F1-model_v4_hydrogens.cif", MMCIFFormat)
specializeresnames!(struc)
chain = struc["A"]
bp, dihedrals = bondparametrization(chain)
plan = jacobianplan(bp)

frame = (chain[1]["N"], chain[1]["CA"], chain[1]["C"])

function objective_and_grad(bp, plan, θ, frame, target)
    X = atomcoordinates(bp, θ, frame)
    w = X .- target                     # gradient of (1/2) Σ |X - target|² w.r.t. X
    val = 0.5 * sum(v -> sum(abs2, v), w)
    g = vjp!(zeros(ndihedrals(plan)), plan, X, w)
    return val, g
end

target = atomcoordinates(bp, dihedrals, frame)  # a structure to match
val, g = objective_and_grad(bp, plan, dihedrals, frame, target)
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
  `atomcoordinates`, and the bond angles and fixed dihedrals stored in `bp`.

- **Coordinate layout.** `atomcoordinates`, `jvp!`, `vjp!`, and
  `weightedhessian!` represent a configuration as a
  `Vector{SVector{3,T}}` with one entry per atom, while
  `coordinatejacobian` returns a dense `3*natoms × ndih` matrix whose rows
  `3(t-1)+1:3t` belong to atom `t`. The two representations share memory
  layout: for `X::Vector{SVector{3,T}}`, `reinterpret(T, X)` is the
  corresponding flat vector of length `3*natoms`, and
  `reinterpret(reshape, T, X)` is the same data as a `3 × natoms` matrix.

  A per-atom cotangent `w::Vector{SVector{3,T}}` flattens the same way, so

  ```julia
  J = coordinatejacobian(plan, X)
  J' * reinterpret(T, w) ≈ vjp(plan, X, w)
  J * v                  ≈ reinterpret(T, jvp(plan, X, v))
  ```

  `J` is an ordinary matrix, so `J' * J`, `J \ r`, and factorizations work
  directly.

- **Layout and ordering of `dihedrals`.** `bp.steps` is a single build
  sequence: atoms `1:3` are residue 1's N, Cα, and C, and every later atom is
  placed by exactly one step, in the order the atoms appear in `bp.atoms`.
  A step either extends the chain `a–b–c–d` by one atom `d`, using a bond
  length, a bond angle, and a dihedral about the `b–c` bond, or places
  several atoms at once in a rigid local frame (the internal `Extend` and
  `Branch` types). Only the first kind carries a `rotatable` field, and the
  entries of `dihedrals` (of length `ndihedrals(bp)`) correspond, in order,
  to the steps for which it is `true`.
  [`dihedralangles`](@ref) measures that vector from a set of coordinates,
  inverting `atomcoordinates`. For a chain of `nres` residues the steps run:

  1. **Backbone**, residue by residue from residue 2 to residue `nres`: the
     step placing N_i (rotatable, ψ_{i-1}, about the Cα_{i-1}–C_{i-1} bond),
     the step placing Cα_i (fixed by ω_{i-1}), and the step placing C_i
     (φ_i, about the N_i–Cα_i bond, rotatable unless a ring fixes it). After
     the last residue's C comes the step placing the terminal OXT atom
     (rotatable, about the Cα_nres–C_nres bond).
  2. **Sidechains**, residue by residue from residue 1 to residue `nres`: for
     each residue, its build sequence's own steps in order (typically χ1,
     χ2, … outward from Cα, though a residue's sequence can also expose a
     single-bond rotation that is not conventionally numbered as a χ angle,
     such as a terminal methyl's orientation).

  [`dihedrallabels`](@ref) returns this ordering as data: one
  [`DihedralLabel`](@ref) per entry of `dihedrals`, and hence per column of
  `J`. Each label carries the residue number of the rotation axis, the name
  `:ψ` or `:φ` for a backbone dihedral (`nothing` otherwise), and the four
  `AtomKey` values `a`, `b`, `c`, `d` of the dihedral `a–b–c–d`. Side-chain
  dihedrals are unnamed because the build tables define them from reference
  atoms that need not be the IUPAC χ reference atoms — lysine's rotation
  about Cα–Cβ is measured here as C–Cα–Cβ–Cγ rather than N–Cα–Cβ–Cγ — so the
  four atoms, not a χ number, are what pin down the angle.

- **Fixed vs. rotatable dihedrals.** A step whose `rotatable` field is
  `false`, and every step that places several atoms at once, uses geometry
  stored in `bp` rather than an entry of `dihedrals`. Three kinds of dihedral
  angle are fixed this way:

  + Every ω (the dihedral about each peptide C–N bond), held in the `ϕ`
    field of the step that places Cα.
  + φ for a residue whose ring fixes it, such as proline: the step placing
    that residue's C is marked non-rotatable and holds the fixed value.
    Proline's ring also fixes its χ1 and χ2, whose steps are likewise
    non-rotatable rather than entries in `dihedrals`.
  + Any other non-rotatable build step, whose fixed value is stored in that
    step's own `ϕ` field.
