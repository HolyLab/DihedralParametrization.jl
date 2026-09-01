# Analytic derivatives of `atomcoordinates` with respect to rotatable dihedrals.

# Each rotatable dihedral rotates its descendants about its frame axis:
#
#   ∂X[t]/∂θ_k = u_k × (X[t] - p_k),   u_k = (X[c_k]-X[b_k])/‖X[c_k]-X[b_k]‖,  p_k = X[c_k]

"""
    plan::JacobianPlan

A reusable plan for the coordinate Jacobian ∂X/∂θ and related products.
Construct one with `jacobianplan`.

[`natoms`](@ref) and [`ndihedrals`](@ref) report the sizes of the
coordinate and dihedral vectors the plan expects.

# Extended help

The fields `natoms::Int` and `ndih::Int` back those accessors; the
remaining fields encode the internal build-order tree.
"""
struct JacobianPlan
    natoms::Int
    ndih::Int
    deepest::Vector{Int}
    parent::Vector{Int}
    bidx::Vector{Int}
    cidx::Vector{Int}
end

Base.show(io::IO, plan::JacobianPlan) = print(io, "JacobianPlan with $(plan.natoms) atoms and $(plan.ndih) dihedrals")

Base.:(==)(x::JacobianPlan, y::JacobianPlan) = _fieldsmatch(==, x, y)
Base.isequal(x::JacobianPlan, y::JacobianPlan) = _fieldsmatch(isequal, x, y)
Base.hash(x::JacobianPlan, h::UInt) = _hashfields(x, hash(:JacobianPlan, h))

ndihedrals(plan::JacobianPlan) = plan.ndih
natoms(plan::JacobianPlan) = plan.natoms

# Name one atom of `bp` for an error message.
_atomdesc(bp::BondParametrization, t::Int) =
    (a = bp.atoms[t]; "residue $(a.resnum) $(a.aname)")

# Describe a build step for an error message.
_stepdesc(bp::BondParametrization, aidx::Int, predecessors::Tuple{Int,Int,Int}) =
    "placing $(_atomdesc(bp, aidx)) from " *
    join((_atomdesc(bp, p) for p in predecessors), ", ")

# Return the deepest dihedral moving the frame, checking that its movesets
# are nested and that omitted dihedrals rotate about a frame atom. The
# description of the offending step is built only when an error is raised.
function _framedeepest(deepest::Vector{Int}, parent::Vector{Int}, bidx::Vector{Int}, cidx::Vector{Int},
                       a::Int, b::Int, c::Int, bp::BondParametrization, aidx::Int)
    da, db, dc = deepest[a], deepest[b], deepest[c]
    dmax = max(da, db, dc)
    dmax == 0 && return 0
    for (p, dp) in ((a, da), (b, db), (c, dc))
        k = dmax
        while k != dp
            k == 0 && throw(ArgumentError("$(_stepdesc(bp, aidx, (a, b, c))): coordinate map has a mixed frame the analytic Jacobian cannot represent " *
                                          "(atom $p's moveset is not nested with dihedral $dmax's chain)"))
            if p != bidx[k] && p != cidx[k]
                throw(ArgumentError("$(_stepdesc(bp, aidx, (a, b, c))): coordinate map has a mixed frame the analytic Jacobian cannot represent " *
                                    "(atom $p lacks dihedral $k in its moveset but is not one of $k's axis atoms)"))
            end
            k = parent[k]
        end
    end
    return dmax
end

"""
    plan = jacobianplan(bp::BondParametrization)

Build a derivative plan for coordinates produced from `bp`.

Throws an `ArgumentError` when a build step's frame is incompatible with the
rigid-rotation formula, or when `bp` is internally inconsistent.
"""
function jacobianplan(bp::BondParametrization)
    natoms = length(bp.atoms)
    bp.nres >= 1 || throw(ArgumentError("bp must contain at least one residue"))
    deepest = zeros(Int, natoms)
    parent = Int[]
    bidx = Int[]
    cidx = Int[]

    # Atoms 1:3 are the fixed N, Cα, C reference frame; the rest follow
    # `bp.steps`, exactly as in `atomcoordinates`.
    nplaced = 3
    for step in bp.steps
        nplaced + 1 == step.aidx ||
            throw(ArgumentError("build step places atom $(step.aidx) but $nplaced atoms have been placed; bp.steps is out of order"))
        a, b, c = step.predecessors
        dmax = _framedeepest(deepest, parent, bidx, cidx, a, b, c, bp, step.aidx)
        if step isa Extend
            nplaced += 1
            if step.rotatable
                push!(parent, dmax)
                push!(bidx, b)
                push!(cidx, c)
                deepest[step.aidx] = length(parent)
            else
                deepest[step.aidx] = dmax
            end
        else
            for j in eachindex(step.βs)
                nplaced += 1
                deepest[step.aidx + j - firstindex(step.βs)] = dmax
            end
        end
    end
    nplaced == natoms ||
        throw(ArgumentError("bp.steps placed $nplaced atoms but bp.atoms has $natoms"))

    return JacobianPlan(natoms, ndihedrals(bp), deepest, parent, bidx, cidx)
end

function _axis(plan::JacobianPlan, X::AbstractVector, k::Int)
    b, c = X[plan.bidx[k]], X[plan.cidx[k]]
    bc = c - b
    return bc / norm(bc), c
end

function _checklengths(plan::JacobianPlan, X, name)
    length(X) == plan.natoms ||
        throw(DimensionMismatch("length($name) = $(length(X)) does not match plan's $(plan.natoms) atoms"))
    return nothing
end

"""
    J = coordinatejacobian(plan::JacobianPlan, X::AbstractVector{<:SVector{3}})

Compute ∂X/∂θ as a dense `3 * length(X) × plan.ndih` matrix. The rows for
atom `t` are `3(t-1)+1:3t`. `X` must have been produced from the same bond
parametrization as `plan`.

A coordinate list `Y::Vector{SVector{3,T}}` corresponds to the flat vector
`reinterpret(T, Y)`, so `J * v` matches `reinterpret(T, jvp(plan, X, v))`
and `J' * reinterpret(T, w)` matches `vjp(plan, X, w)`.
"""
function coordinatejacobian(plan::JacobianPlan, X::AbstractVector{<:SVector{3}})
    T = eltype(eltype(X))
    J = zeros(T, 3 * plan.natoms, plan.ndih)
    return coordinatejacobian!(J, plan, X)
end

"""
    coordinatejacobian!(J, plan::JacobianPlan, X::AbstractVector{<:SVector{3}})

In-place version of `coordinatejacobian`. `J` is overwritten and must have
size `(3 * length(X), plan.ndih)`.
"""
function coordinatejacobian!(J::AbstractMatrix, plan::JacobianPlan, X::AbstractVector{<:SVector{3}})
    _checklengths(plan, X, "X")
    size(J) == (3 * plan.natoms, plan.ndih) ||
        throw(DimensionMismatch("size(J) = $(size(J)) does not match plan's $((3 * plan.natoms, plan.ndih))"))
    fill!(J, zero(eltype(J)))
    for t in eachindex(X)
        k = plan.deepest[t]
        k == 0 && continue
        xt = X[t]
        rows = 3 * (t - 1) + 1 : 3 * t
        while k != 0
            u, p = _axis(plan, X, k)
            J[rows, k] = cross(u, xt - p)
            k = plan.parent[k]
        end
    end
    return J
end

"""
    g = vjp!(g, plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, w::AbstractVector{<:SVector{3}})

Compute `g = J' * w` without forming `J`, where `J = ∂X/∂θ`. `w[t]` is the
weight for atom `t`. Returns `g`.
"""
function vjp!(g::AbstractVector, plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, w::AbstractVector{<:SVector{3}})
    _checklengths(plan, X, "X")
    _checklengths(plan, w, "w")
    length(g) == plan.ndih || throw(DimensionMismatch("length(g) = $(length(g)) does not match plan's $(plan.ndih) dihedrals"))
    T = promote_type(eltype(eltype(X)), eltype(eltype(w)))
    S0 = fill(zero(SVector{3,T}), plan.ndih)
    S1 = fill(zero(SVector{3,T}), plan.ndih)
    for t in eachindex(X)
        k = plan.deepest[t]
        k == 0 && continue
        S0[k] += w[t]
        S1[k] += cross(X[t], w[t])
    end
    # Parents precede children in build order.
    for k = plan.ndih:-1:1
        pk = plan.parent[k]
        if pk != 0
            S0[pk] += S0[k]
            S1[pk] += S1[k]
        end
    end
    for k = 1:plan.ndih
        u, p = _axis(plan, X, k)
        g[k] = dot(u, S1[k] - cross(p, S0[k]))
    end
    return g
end

"""
    g = vjp(plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, w::AbstractVector{<:SVector{3}})

Allocating version of `vjp!`.
"""
function vjp(plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, w::AbstractVector{<:SVector{3}})
    T = promote_type(eltype(eltype(X)), eltype(eltype(w)))
    g = Vector{T}(undef, plan.ndih)
    return vjp!(g, plan, X, w)
end

"""
    δx = jvp!(δx, plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, v::AbstractVector)

Compute `δx = J * v` without forming `J`, where `J = ∂X/∂θ`. `v` has one
entry per dihedral and `δx` one entry per atom. Returns `δx`.
"""
function jvp!(δx::AbstractVector, plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, v::AbstractVector)
    _checklengths(plan, X, "X")
    length(δx) == plan.natoms || throw(DimensionMismatch("length(δx) = $(length(δx)) does not match plan's $(plan.natoms) atoms"))
    length(v) == plan.ndih || throw(DimensionMismatch("length(v) = $(length(v)) does not match plan's $(plan.ndih) dihedrals"))
    T = promote_type(eltype(eltype(X)), eltype(v))
    Ω = Vector{SVector{3,T}}(undef, plan.ndih)
    τ = Vector{SVector{3,T}}(undef, plan.ndih)
    for k = 1:plan.ndih
        u, p = _axis(plan, X, k)
        pk = plan.parent[k]
        Ωpar = pk == 0 ? zero(SVector{3,T}) : Ω[pk]
        τpar = pk == 0 ? zero(SVector{3,T}) : τ[pk]
        Ω[k] = Ωpar + v[k] * u
        τ[k] = τpar - v[k] * cross(u, p)
    end
    for t in eachindex(X)
        k = plan.deepest[t]
        δx[t] = k == 0 ? zero(SVector{3,T}) : cross(Ω[k], X[t]) + τ[k]
    end
    return δx
end

"""
    δx = jvp(plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, v::AbstractVector)

Allocating version of `jvp!`.
"""
function jvp(plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, v::AbstractVector)
    T = promote_type(eltype(eltype(X)), eltype(v))
    δx = Vector{SVector{3,T}}(undef, plan.natoms)
    return jvp!(δx, plan, X, v)
end

"""
    S = weightedhessian!(S, plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, w::AbstractVector{<:SVector{3}})

Compute `S[i,j] = Σ_a w[a] ⋅ ∂²X[a]/∂θ_i∂θ_j`, the Hessian with respect to
the dihedrals θ of the scalar `Σ_a w[a] ⋅ X[a](θ)`, evaluated at the
configuration `X`. This is not a Hessian-vector product: `S` is the full
`ndih × ndih` Hessian of `w ⋅ X`, computed without forming the coordinate
map's second-derivative tensor.

`X` and `w` each have one `SVector{3}` entry per atom, `w[t]` playing the
same role as in `vjp!`. `S` is symmetric and is overwritten; it must have
size `(plan.ndih, plan.ndih)`. Returns `S`.
"""
function weightedhessian!(S::AbstractMatrix, plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, w::AbstractVector{<:SVector{3}})
    _checklengths(plan, X, "X")
    _checklengths(plan, w, "w")
    size(S) == (plan.ndih, plan.ndih) ||
        throw(DimensionMismatch("size(S) = $(size(S)) does not match plan's $((plan.ndih, plan.ndih)) dihedrals"))
    fill!(S, zero(eltype(S)))
    T = promote_type(eltype(eltype(X)), eltype(eltype(w)))
    S0 = fill(zero(SVector{3,T}), plan.ndih)
    s = zeros(T, plan.ndih)
    M = fill(zero(SMatrix{3,3,T,9}), plan.ndih)
    for t in eachindex(X)
        k = plan.deepest[t]
        k == 0 && continue
        S0[k] += w[t]
        s[k] += dot(X[t], w[t])
        M[k] += X[t] * w[t]'
    end
    # Parents precede children in build order.
    for k = plan.ndih:-1:1
        pk = plan.parent[k]
        if pk != 0
            S0[pk] += S0[k]
            s[pk] += s[k]
            M[pk] += M[k]
        end
    end
    u = Vector{SVector{3,T}}(undef, plan.ndih)
    p = Vector{SVector{3,T}}(undef, plan.ndih)
    for k = 1:plan.ndih
        u[k], p[k] = _axis(plan, X, k)
    end
    for j = 1:plan.ndih
        uj, pj = u[j], p[j]
        Tj = M[j] * uj - pj * dot(uj, S0[j]) - uj * (s[j] - dot(pj, S0[j]))
        i = j
        while i != 0
            Sij = dot(u[i], Tj)
            S[i, j] = Sij
            S[j, i] = Sij
            i = plan.parent[i]
        end
    end
    return S
end

"""
    S = weightedhessian(plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, w::AbstractVector{<:SVector{3}})

Allocating version of `weightedhessian!`.
"""
function weightedhessian(plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, w::AbstractVector{<:SVector{3}})
    T = promote_type(eltype(eltype(X)), eltype(eltype(w)))
    S = zeros(T, plan.ndih, plan.ndih)
    return weightedhessian!(S, plan, X, w)
end
