## Analytic Jacobian of `atomcoordinates` with respect to the rotatable dihedrals.

# Each rotatable dihedral rotates its descendants about its frame axis:
#
#   ∂X[t]/∂θ_k = u_k × (X[t] - p_k),   u_k = (X[c_k]-X[b_k])/‖X[c_k]-X[b_k]‖,  p_k = X[c_k]

"""
    plan::JacobianPlan

A reusable plan for evaluating the coordinate-map Jacobian. Construct one
with `jacobianplan`.
"""
struct JacobianPlan
    natoms::Int
    ndih::Int
    deepest::Vector{Int}
    parent::Vector{Int}
    bidx::Vector{Int}
    cidx::Vector{Int}
end

# Return the deepest dihedral moving the frame, checking that its movesets
# are nested and that omitted dihedrals rotate about a frame atom.
function _framedeepest(deepest::Vector{Int}, parent::Vector{Int}, bidx::Vector{Int}, cidx::Vector{Int},
                        a::Int, b::Int, c::Int, desc::AbstractString)
    da, db, dc = deepest[a], deepest[b], deepest[c]
    dmax = max(da, db, dc)
    dmax == 0 && return 0
    for (p, dp) in ((a, da), (b, db), (c, dc))
        k = dmax
        while k != dp
            k == 0 && error("$desc: coordinate map has a mixed frame the analytic Jacobian cannot represent " *
                             "(atom $p's moveset is not nested with dihedral $dmax's chain)")
            if p != bidx[k] && p != cidx[k]
                error("$desc: coordinate map has a mixed frame the analytic Jacobian cannot represent " *
                      "(atom $p lacks dihedral $k in its moveset but is not one of $k's axis atoms)")
            end
            k = parent[k]
        end
    end
    return dmax
end

"""
    plan = jacobianplan(bp::BondParametrization)

Build a reusable plan for evaluating the coordinate-map Jacobian. The plan
may be used with any coordinates produced from the same `bp`.

Throws an error when a build step's frame is incompatible with the
rigid-rotation formula.
"""
function jacobianplan(bp::BondParametrization)
    natoms = length(bp.atoms)
    nres = length(bp.residues)
    deepest = zeros(Int, natoms)
    parent = Int[]
    bidx = Int[]
    cidx = Int[]

    function newdihedral!(b::Int, c::Int, dmax::Int)
        push!(parent, dmax)
        push!(bidx, b)
        push!(cidx, c)
        return length(parent)
    end

    # Fixed N, Cα, C reference frame.
    aidx = 3
    nres >= 1 || error("bp must contain at least one residue")
    prev3, prev2, prev1 = 1, 2, 3
    for i = 2:nres
        # N_i: ψ_{i-1} about Cα_{i-1}–C_{i-1}.
        dmax = _framedeepest(deepest, parent, bidx, cidx, prev3, prev2, prev1, "residue $i N (ψ_$(i-1))")
        aidx += 1
        k = newdihedral!(prev2, prev1, dmax)
        deepest[aidx] = k
        prev3, prev2, prev1 = prev2, prev1, aidx
        # Cα_i: fixed ω_{i-1}.
        dmax = _framedeepest(deepest, parent, bidx, cidx, prev3, prev2, prev1, "residue $i Cα (ω)")
        aidx += 1
        deepest[aidx] = dmax
        prev3, prev2, prev1 = prev2, prev1, aidx
        # C_i: φ_i about N_i–Cα_i when rotatable.
        dmax = _framedeepest(deepest, parent, bidx, cidx, prev3, prev2, prev1, "residue $i C (φ_$i)")
        aidx += 1
        if bp.phirotatable[i]
            k = newdihedral!(prev2, prev1, dmax)
            deepest[aidx] = k
        else
            deepest[aidx] = dmax
        end
        prev3, prev2, prev1 = prev2, prev1, aidx
    end
    # OXT: final ψ′ about Cα_nres–C_nres.
    dmax = _framedeepest(deepest, parent, bidx, cidx, prev3, prev2, prev1, "OXT (ψ′)")
    aidx += 1
    k = newdihedral!(prev2, prev1, dmax)
    deepest[aidx] = k

    for (i, r) in enumerate(bp.residues)
        for step in r.steps
            a, b, c = step.predecessors
            if step isa Extend
                dmax = _framedeepest(deepest, parent, bidx, cidx, a, b, c,
                                      "residue $i sidechain Extend $(step.predecessors)")
                aidx += 1
                if step.rotatable
                    k = newdihedral!(b, c, dmax)
                    deepest[aidx] = k
                else
                    deepest[aidx] = dmax
                end
            elseif step isa Branch
                dmax = _framedeepest(deepest, parent, bidx, cidx, a, b, c,
                                      "residue $i sidechain Branch $(step.predecessors)")
                for _ = 1:length(step.βs)
                    aidx += 1
                    deepest[aidx] = dmax
                end
            else
                error("residue $i: unrecognized build step $step")
            end
        end
    end
    aidx == natoms || error("built $aidx atoms but bp.atoms has $natoms; traversal does not match atomcoordinates")

    ndih = length(parent)
    nphi_free = count(view(bp.phirotatable, 2:nres))
    nsidechain_free = sum((count(step -> step isa Extend && step.rotatable, r.steps) for r in bp.residues); init=0)
    expected = (nres - 1) + nphi_free + 1 + nsidechain_free
    ndih == expected || error("plan has $ndih dihedral columns but atomcoordinates would consume $expected; " *
                               "traversal does not match atomcoordinates")

    return JacobianPlan(natoms, ndih, deepest, parent, bidx, cidx)
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
    g = jtv!(g, plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, w::AbstractVector{<:SVector{3}})

Compute `g = J' * w` without forming `J`, where `J = ∂X/∂θ`. `w[t]` is the
weight for atom `t`. Returns `g`.
"""
function jtv!(g::AbstractVector, plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, w::AbstractVector{<:SVector{3}})
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
    S = vhp!(S, plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, w::AbstractVector{<:SVector{3}})

Compute `S[i,j] = Σ_a w[a] ⋅ ∂²X[a]/∂θ_i∂θ_j`, the contraction of the
coordinate map's second derivative with a per-atom cotangent `w`, without
forming the second-derivative tensor. `X` and `w` each have one `SVector{3}`
entry per atom, `w[t]` playing the same role as in `jtv!`; both the second
derivative and `S` are evaluated at the fixed configuration `X`. `S` is
symmetric and is overwritten; it must have size `(plan.ndih, plan.ndih)`.
Returns `S`.
"""
function vhp!(S::AbstractMatrix, plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, w::AbstractVector{<:SVector{3}})
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
    S = vhp(plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, w::AbstractVector{<:SVector{3}})

Allocating version of `vhp!`.
"""
function vhp(plan::JacobianPlan, X::AbstractVector{<:SVector{3}}, w::AbstractVector{<:SVector{3}})
    T = promote_type(eltype(eltype(X)), eltype(eltype(w)))
    S = zeros(T, plan.ndih, plan.ndih)
    return vhp!(S, plan, X, w)
end
