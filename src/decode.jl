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

function add_to_middle!(X, aidx::Int, a::AbstractVector, b::AbstractVector, c::AbstractVector, βs)
    # Place atoms such as Cβ in tetrahedral geometry, at `X[aidx:aidx+length(βs)-1]`.
    ab = b - a
    ab = ab / norm(ab)
    cb = b - c
    cb = cb / norm(cb)
    n = cross(ab, cb)
    n = n / norm(n)
    M = [ab cb n]
    for (i, β) in pairs(βs)
        X[aidx + i - firstindex(βs)] = b + M * β
    end
    return X
end

add_to_middle(a::AbstractVector, b::AbstractVector, c::AbstractVector, βs) =
    add_to_middle!(Vector{promote_type(typeof(a), typeof(b), typeof(c))}(undef, length(βs)), 1, a, b, c, βs)

# Convert a reference frame to a homogeneous tuple of static 3-vectors.
frametuple(frame::NTuple{3,SVector{3,T}}) where {T<:Real} = frame
function frametuple((n, cα, c)::NTuple{3,AbstractVector})
    T = promote_type(eltype(n), eltype(cα), eltype(c))
    return (SVector{3,T}(n), SVector{3,T}(cα), SVector{3,T}(c))
end
frametuple(frame::NTuple{3,Atom}) = frametuple(map(a -> a.coords, frame))
function frametuple(chain::Chain)
    nter = first(chain)::Residue
    return frametuple((nter["N"]::Atom, nter["CA"]::Atom, nter["C"]::Atom))
end

# The frame with residue 1's N at the origin, Cα on the positive x-axis, and
# C in the xy-plane with positive y.
function canonicalframe(bp::BondParametrization{T}) where {T}
    sθ, cθ = sincos(bp.θncac)
    n = zero(SVector{3,T})
    cα = SVector{3,T}(bp.ℓnca, zero(T), zero(T))
    c = SVector{3,T}(bp.ℓnca - bp.ℓcac * cθ, bp.ℓcac * sθ, zero(T))
    return (n, cα, c)
end

"""
    X = atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector, (n, cα, c))
    X = atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector, chain::Chain)
    X = atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector)

Compute atom coordinates from `bp`, rotatable angles `dihedrals` (in radians),
and a reference frame given by residue 1's backbone N, Cα, and C. Results
follow the order of `bp.atoms`.

The frame may be a tuple `(n, cα, c)` of three 3-vectors or three `Atom`s,
or a `Chain` whose first residue supplies the atoms. The two-argument form
places residue 1's N at the origin, Cα on the positive x-axis, and C in the
xy-plane with positive y.

`length(dihedrals)` must equal `ndihedrals(bp)`; a `DimensionMismatch` is
thrown otherwise. Reference coordinates whose N–Cα distance, Cα–C distance,
or N–Cα–C angle disagrees with `bp` raise an `ArgumentError`, as does a `bp`
whose build sequence is inconsistent with its own atom count.
`X` is a `Vector{SVector{3,R}}`, where `R` promotes the element types of
`bp`, `dihedrals`, and the reference coordinates. This permits automatic
differentiation types in `dihedrals`.

See [`atomcoordinates!`](@ref) for the in-place form.

# Extended help

Reference coordinates are converted to `SVector{3,T}` using their promoted
element type. A tuple of `Atom`s contributes each atom's `coords`; the
`Chain` method uses the first residue's N, CA, and C atoms.
"""
function atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector, frame)
    n, cα, c = frametuple(frame)
    R = promote_type(eltype(bp), eltype(dihedrals), eltype(n))
    X = Vector{SVector{3,R}}(undef, length(bp.atoms))
    return _atomcoordinates!(X, bp, dihedrals, n, cα, c)
end
atomcoordinates(bp::BondParametrization, dihedrals::AbstractVector) =
    atomcoordinates(bp, dihedrals, canonicalframe(bp))

"""
    atomcoordinates!(X, bp::BondParametrization, dihedrals::AbstractVector, (n, cα, c))
    atomcoordinates!(X, bp::BondParametrization, dihedrals::AbstractVector, chain::Chain)
    atomcoordinates!(X, bp::BondParametrization, dihedrals::AbstractVector)

Write [`atomcoordinates`](@ref) into `X` and return it. `X` must contain
`length(bp.atoms)` 3-vectors. Coordinates are converted to `eltype(X)`.
"""
atomcoordinates!(X::AbstractVector, bp::BondParametrization, dihedrals::AbstractVector, frame) =
    _atomcoordinates!(X, bp, dihedrals, frametuple(frame)...)
atomcoordinates!(X::AbstractVector, bp::BondParametrization, dihedrals::AbstractVector) =
    atomcoordinates!(X, bp, dihedrals, canonicalframe(bp))

function _atomcoordinates!(X::AbstractVector, bp::BondParametrization, dihedrals::AbstractVector,
                           n::SVector{3,Tref}, cα::SVector{3,Tref}, c::SVector{3,Tref}) where {Tref<:Real}
    Base.require_one_based_indexing(X)
    # Check that the inputs are consistent with `bp`
    nd = ndihedrals(bp)
    length(dihedrals) == nd ||
        throw(DimensionMismatch("length(dihedrals) = $(length(dihedrals)) does not match bp's $nd rotatable dihedrals"))
    natoms = length(bp.atoms)
    length(X) == natoms ||
        throw(DimensionMismatch("length(X) = $(length(X)) does not match bp's $natoms atoms"))
    ℓnca, ℓcac = norm(cα - n), norm(c - cα)
    θncac = bondangle(n - cα, c - cα)
    ℓnca ≈ bp.ℓnca ||
        throw(ArgumentError("the reference N–Cα distance is $ℓnca but bp requires $(bp.ℓnca)"))
    ℓcac ≈ bp.ℓcac ||
        throw(ArgumentError("the reference Cα–C distance is $ℓcac but bp requires $(bp.ℓcac)"))
    θncac ≈ bp.θncac ||
        throw(ArgumentError("the reference N–Cα–C angle is $θncac but bp requires $(bp.θncac)"))

    # Atoms 1:3 are the reference frame; `bp.steps` places all the rest.
    R = promote_type(eltype(bp), eltype(dihedrals), Tref)
    X[1], X[2], X[3] = n, cα, c
    placed = 3                        # number of atoms placed so far
    idx = firstindex(dihedrals) - 1   # index of the last dihedral consumed
    for step in bp.steps
        placed + 1 == step.aidx ||
            throw(ArgumentError("build step places atom $(step.aidx) but $placed atoms have been placed; bp.steps is out of order"))
        nplaced = step isa Extend ? 1 : length(step.βs)
        placed + nplaced <= natoms ||
            throw(ArgumentError("bp.steps places atom $(placed + nplaced) but bp.atoms has $natoms"))
        a, b, cc = X[SVector(step.predecessors)]
        if step isa Extend
            ϕ = step.rotatable ? convert(R, dihedrals[idx+=1]) : convert(R, step.ϕ)
            X[step.aidx] = snnerf(a, b, cc, step.ℓcd, step.θbcd, ϕ)
        else
            add_to_middle!(X, step.aidx, a, b, cc, step.βs)
        end
        placed += nplaced
    end
    placed == natoms ||
        throw(ArgumentError("bp.steps placed $placed atoms but bp.atoms has $natoms"))
    return X
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
buildchain(reference::Chain, bp::BondParametrization, X::AbstractVector{<:AbstractVector}) =
    buildchain(reference, bp, svectors(X))
