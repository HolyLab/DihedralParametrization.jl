"""
    d = snnerf(a, b, c, ℓcd, θbcd, ϕbc)

Compute point `d` from three preceding points and its bond length, bond angle,
and dihedral angle using the SN-NeRF algorithm. Angles are in radians.

## Reference

> Parsons, Jerod, et al. "Practical conversion from torsion space to Cartesian
> space for in silico protein synthesis." Journal of computational chemistry 26.10 (2005): 1063-1068.
"""
function snnerf(a::AbstractVector, b::AbstractVector, c::AbstractVector,
                ℓcd::Real, θbcd::Real, ϕbc::Real)
    return snnerf(a, b, c, ℓcd, sincos(θbcd), sincos(ϕbc))
end

function snnerf(a::AbstractVector, b::AbstractVector, c::AbstractVector,
                ℓcd::Real, θbcd::Real, (sϕ, cϕ)::Tuple{Real,Real})
    return snnerf(a, b, c, ℓcd, sincos(θbcd), (sϕ, cϕ))
end

function snnerf(a::AbstractVector, b::AbstractVector, c::AbstractVector,
                ℓcd::Real, (sθ, cθ)::Tuple{Real,Real}, (sϕ, cϕ)::Tuple{Real,Real})
    bc = c - b
    bc = bc / norm(bc)
    ab = b - a
    n = cross(ab, bc)
    n = n / norm(n)
    M = [bc cross(n, bc) n]
    ℓsθ, ℓcθ = ℓcd .* (sθ, cθ)
    d2 = SVector(-ℓcθ, ℓsθ * cϕ, ℓsθ * sϕ)
    return c + M * d2
end

function add_to_middle!(X, a::AbstractVector, b::AbstractVector, c::AbstractVector, βs)
    # Place atoms such as Cβ in tetrahedral geometry.
    ab = b - a
    ab = ab / norm(ab)
    cb = b - c
    cb = cb / norm(cb)
    n = cross(ab, cb)
    n = n / norm(n)
    M = [ab cb n]
    for β in βs
        push!(X, b + M * β)
    end
    return X
end

add_to_middle(a::AbstractVector, b::AbstractVector, c::AbstractVector, βs) =
    add_to_middle!(promote_type(typeof(a), typeof(b), typeof(c))[], a, b, c, βs)

"""
    X = atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector, n::AbstractVector, cα::AbstractVector, c::AbstractVector)
    X = atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector, n::Atom, cα::Atom, c::Atom)
    X = atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector, chain::Chain)
    X = atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector)

Compute atom coordinates from `bp`, rotatable angles `dihedrals` (in radians),
and the first backbone N, Cα, and C coordinates. Results follow the order of
`bp.atoms`.

The two-argument form places residue 1's N at the origin, Cα on the positive
x-axis, and C in the xy-plane with positive y.

`length(dihedrals)` must equal `ndihedrals(bp)`; a `DimensionMismatch` is
thrown otherwise. Reference coordinates whose N–Cα distance, Cα–C distance,
or N–Cα–C angle disagrees with `bp` raise an `ArgumentError`, as does a `bp`
whose build sequence is inconsistent with its own atom and dihedral counts.
`X` is a `Vector{SVector{3,R}}`, where `R` promotes the element types of
`bp`, `dihedrals`, and the reference coordinates. This permits automatic
differentiation types in `dihedrals`.

# Extended help

Reference coordinates are converted to `SVector{3,T}` using their promoted
element type. The `Atom` method reads `coords`; the `Chain` method uses the
first residue's N, CA, and C atoms.
"""
function atomcoordinates(bp::BondParametrization{T}, dihedrals::AbstractVector{S},
                         n::SVector{3,Tref}, cα::SVector{3,Tref}, c::SVector{3,Tref}) where {T<:Real, S<:Real, Tref<:Real}
    # Check that the inputs are consistent with `bp`
    nd = ndihedrals(bp)
    length(dihedrals) == nd ||
        throw(DimensionMismatch("length(dihedrals) = $(length(dihedrals)) does not match bp's $nd rotatable dihedrals"))
    ℓnca, ℓcac = norm(cα - n), norm(c - cα)
    θncac = bondangle(n - cα, c - cα)
    ℓnca ≈ bp.ℓnca ||
        throw(ArgumentError("the reference N–Cα distance is $ℓnca but bp requires $(bp.ℓnca)"))
    ℓcac ≈ bp.ℓcac ||
        throw(ArgumentError("the reference Cα–C distance is $ℓcac but bp requires $(bp.ℓcac)"))
    θncac ≈ bp.θncac ||
        throw(ArgumentError("the reference N–Cα–C angle is $θncac but bp requires $(bp.θncac)"))

    # Atoms 1:3 are the reference frame; `bp.steps` places all the rest.
    R = promote_type(T, S, Tref)
    X = sizehint!(SVector{3,R}[n, cα, c], length(bp.atoms))
    idx = firstindex(dihedrals) - 1   # index of the last dihedral consumed
    for step in bp.steps
        length(X) + 1 == step.aidx ||
            throw(ArgumentError("build step places atom $(step.aidx) but $(length(X)) atoms have been placed; bp.steps is out of order"))
        a, b, cc = X[SVector(step.predecessors)]
        if step isa Extend
            ϕ = step.rotatable ? convert(R, dihedrals[idx+=1]) : convert(R, step.ϕ)
            push!(X, snnerf(a, b, cc, step.ℓcd, step.θbcd, ϕ))
        else
            add_to_middle!(X, a, b, cc, step.βs)
        end
    end
    length(X) == length(bp.atoms) ||
        throw(ArgumentError("bp.steps placed $(length(X)) atoms but bp.atoms has $(length(bp.atoms))"))
    nconsumed = idx - firstindex(dihedrals) + 1
    nconsumed == nd ||
        throw(ArgumentError("bp.steps consumed $nconsumed dihedrals but bp declares $nd"))
    return X
end
function atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector, n::AbstractVector, cα::AbstractVector, c::AbstractVector)
    T = promote_type(eltype(n), eltype(cα), eltype(c))
    return atomcoordinates(bp, dihedrals, SVector{3,T}(n), SVector{3,T}(cα), SVector{3,T}(c))
end
atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector, n::Atom, cα::Atom, c::Atom) = atomcoordinates(bp, dihedrals, n.coords, cα.coords, c.coords)
function atomcoordinates(bp::BondParametrization{T}, dihedrals::AbstractVector) where {T<:Real}
    sθ, cθ = sincos(bp.θncac)
    n = zero(SVector{3,T})
    cα = SVector{3,T}(bp.ℓnca, zero(T), zero(T))
    c = SVector{3,T}(bp.ℓnca - bp.ℓcac * cθ, bp.ℓcac * sθ, zero(T))
    return atomcoordinates(bp, dihedrals, n, cα, c)
end
function atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector, chain::Chain)
    nter = first(chain)::Residue
    return atomcoordinates(bp, dihedrals, nter["N"]::Atom, nter["CA"]::Atom, nter["C"]::Atom)
end

"""
    out = buildchain(reference::Chain, bp::BondParametrization, X::AbstractVector{<:SVector{3}})

Copy `reference` and replace its atom coordinates with `X`, matched through
`bp.atoms`. `reference` comes first because it is the template being copied
and overwritten, like the destination argument of `copyto!`.

`length(X)` must equal `length(bp.atoms)`; a `DimensionMismatch` is thrown
otherwise. An atom of `reference` that `bp.atoms` does not name raises an
`ArgumentError`.
"""
function buildchain(reference::Chain, bp::BondParametrization, X::AbstractVector{<:SVector{3}})
    length(X) == length(bp.atoms) ||
        throw(DimensionMismatch("length(X) = $(length(X)) does not match bp's $(length(bp.atoms)) atoms"))
    out = copy(reference)
    coordidx = Dict{AtomKey, Int}()
    for (i, akey) in pairs(bp.atoms)
        coordidx[akey] = i
    end
    for a in collectatoms(out)
        akey = AtomKey(a)
        i = get(coordidx, akey, nothing)
        i === nothing &&
            throw(ArgumentError("reference has atom $(akey.aname) in residue $(akey.resnum), which the parametrization does not describe"))
        a.coords .= X[i]
    end
    return out
end
