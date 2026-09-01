"""
    d = snnerf(a, b, c, ℓcd, θbcd, ϕbc)

Given three points `a`, `b`, and `c` (as 3D vectors), the bond length `ℓcd`,
the bond angle `θbcd` (in radians), and the dihedral angle `ϕbc` (in radians),
compute the coordinates of point `d` using the SN-NeRF algorithm.

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
    # Calculate unit vectors
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
    # used, e.g., to place the Cβ atom in a tetrahedral geometry
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

Compute atom coordinates from `bp`, rotatable angles `dihedrals` (in radians),
and the first backbone N, Cα, and C coordinates. Results follow the order of
`bp.atoms`.

`length(dihedrals)` must equal `ndihedrals(bp)`; a `DimensionMismatch` is
thrown otherwise. The element type of `dihedrals` is independent of the
element type of `bp` and the reference coordinates, so automatic
differentiation can pass dual or tangent numbers as `dihedrals` while `bp`
stays `Float64`. `X` is a `Vector{SVector{3,R}}` with `R` the
`promote_type` of `bp`'s element type, the element type of `dihedrals`, and
the element type of the reference coordinates.

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
    norm(cα - n) ≈ bp.ℓnca || error("Provided N and Cα do not match bond length in bp")
    norm(c - cα) ≈ bp.ℓcac || error("Provided Cα and C do not match bond length in bp")
    bondangle(n - cα, c - cα) ≈ bp.θncac || error("Provided N, Cα, and C do not match bond angle in bp")

    # Atoms 1:3 are the reference frame; `bp.steps` places all the rest.
    R = promote_type(T, S, Tref)
    X = sizehint!(SVector{3,R}[n, cα, c], length(bp.atoms))
    idx = firstindex(dihedrals) - 1   # index of the last dihedral consumed
    for step in bp.steps
        length(X) + 1 == step.aidx ||
            error("build step places atom $(step.aidx) but $(length(X)) atoms have been placed; bp.steps is out of order")
        a, b, cc = X[SVector(step.predecessors)]
        if step isa Extend
            ϕ = step.rotatable ? convert(R, dihedrals[idx+=1]) : convert(R, step.ϕ)
            push!(X, snnerf(a, b, cc, step.ℓcd, step.θbcd, ϕ))
        else
            add_to_middle!(X, a, b, cc, step.βs)
        end
    end
    length(X) == length(bp.atoms) ||
        error("bp.steps placed $(length(X)) atoms but bp.atoms has $(length(bp.atoms))")
    nconsumed = idx - firstindex(dihedrals) + 1
    nconsumed == nd ||
        error("bp.steps consumed $nconsumed dihedrals but bp declares $nd")
    return X
end
function atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector, n::AbstractVector, cα::AbstractVector, c::AbstractVector)
    T = promote_type(eltype(n), eltype(cα), eltype(c))
    return atomcoordinates(bp, dihedrals, SVector{3,T}(n), SVector{3,T}(cα), SVector{3,T}(c))
end
atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector, n::Atom, cα::Atom, c::Atom) = atomcoordinates(bp, dihedrals, n.coords, cα.coords, c.coords)
function atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector, chain::Chain)
    nter = first(chain)::Residue
    return atomcoordinates(bp, dihedrals, nter["N"]::Atom, nter["CA"]::Atom, nter["C"]::Atom)
end

"""
    out = buildchain(reference::Chain, bp::BondParametrization, X::AbstractVector{<:SVector{3}})

Copy `reference` and replace its atom coordinates with `X`, matched through
`bp.atoms`. `reference` comes first because it is the template being copied
and overwritten, like the destination argument of `copyto!`.
"""
function buildchain(reference::Chain, bp::BondParametrization, X::AbstractVector{<:SVector{3}})
    out = copy(reference)
    coordidx = Dict{AtomKey, Int}()
    for (i, akey) in enumerate(bp.atoms)
        coordidx[akey] = i
    end
    for a in collectatoms(out)
        i = coordidx[AtomKey(a)]
        a.coords .= X[i]
    end
    return out
end
