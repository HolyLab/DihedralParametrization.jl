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
    X = atomcoordinates(bp::BondParametrization, dihedrals::Vector, n::AbstractVector, cα::AbstractVector, c::AbstractVector)
    X = atomcoordinates(bp::BondParametrization, dihedrals::Vector, n::Atom, cα::Atom, c::Atom)
    X = atomcoordinates(bp::BondParametrization, dihedrals::Vector, chain::Chain)

Compute atom coordinates from `bp`, rotatable angles `dihedrals` (in radians),
and the first backbone N, Cα, and C coordinates. Results follow the order of
`bp.atoms`.

`length(dihedrals)` must equal `ndihedrals(bp)`; a `DimensionMismatch` is
thrown otherwise. The element type of `dihedrals` is independent of the
element type of `bp` and the reference coordinates, so automatic
differentiation can pass dual or tangent numbers as `dihedrals` while `bp`
stays `Float64`.

# Extended help

Reference coordinates are converted to `SVector{3,T}` using their promoted
element type. The `Atom` method reads `coords`; the `Chain` method uses the
first residue's N, CA, and C atoms.
"""
function atomcoordinates(bp::BondParametrization, dihedrals::Vector{S}, n::SVector{3,T}, cα::SVector{3,T}, c::SVector{3,T}) where {S<:Real, T<:Real}
    # Check that the inputs are consistent with `bp`
    nd = ndihedrals(bp)
    length(dihedrals) == nd ||
        throw(DimensionMismatch("length(dihedrals) = $(length(dihedrals)) does not match bp's $nd rotatable dihedrals"))
    norm(cα - n) ≈ bp.bblengths[1] || error("Provided N and Cα do not match bond length in bp")
    norm(c - cα) ≈ bp.bblengths[2] || error("Provided Cα and C do not match bond length in bp")
    bondangle(n - cα, c - cα) ≈ bp.bbangles[1] || error("Provided N, Cα, and C do not match bond angle in bp")

    # Initialize the coordinates array
    X = sizehint!(SVector{3,promote_type(S, T)}[n, cα, c], length(bp.atoms))
    idx = 0  # index into dihedrals vector
    # Connect the backbone
    prev3, prev2, prev1 = n, cα, c
    nres = length(bp.residues)
    for i = 2:nres
        # add the next N
        ψ = dihedrals[idx+=1]       # rotatable
        ℓ = bp.bblengths[3*i - 3]
        θ = bp.bbangles[3*i - 4]
        d = snnerf(prev3, prev2, prev1, ℓ, θ, ψ)
        push!(X, d)
        prev3, prev2, prev1 = prev2, prev1, d
        # add the next Cα
        ω = bp.omegas[i - 1]
        ℓ = bp.bblengths[3*i - 2]
        θ = bp.bbangles[3*i - 3]
        d = snnerf(prev3, prev2, prev1, ℓ, θ, ω)
        push!(X, d)
        prev3, prev2, prev1 = prev2, prev1, d
        # add the next C
        φ = bp.phirotatable[i] ? dihedrals[idx+=1] : bp.phi[i]
        ℓ = bp.bblengths[3*i - 1]
        θ = bp.bbangles[3*i - 2]
        d = snnerf(prev3, prev2, prev1, ℓ, θ, φ)
        push!(X, d)
        prev3, prev2, prev1 = prev2, prev1, d
    end
    # add a terminal O
    ψ′ = dihedrals[idx+=1]  # for placement of OXT
    ℓ = bp.bblengths[end]
    θ = bp.bbangles[end]
    d = snnerf(prev3, prev2, prev1, ℓ, θ, ψ′)
    push!(X, d)
    for (i, r) in enumerate(bp.residues)
        for step in r.steps
            if step isa Extend{T}
                a, b, c = X[SVector(step.predecessors)]
                ϕ = step.rotatable ? dihedrals[idx+=1] : step.ϕ
                d = snnerf(a, b, c, step.ℓcd, step.θbcd, ϕ)
                push!(X, d)
            elseif step isa Branch{T}
                a, b, c = X[SVector(step.predecessors)]
                d = add_to_middle!(X, a, b, c, step.βs)
            end
        end
    end
    @assert length(X) == length(bp.atoms)
    return X
end
function atomcoordinates(bp::BondParametrization, dihedrals::Vector, n::AbstractVector, cα::AbstractVector, c::AbstractVector)
    T = promote_type(eltype(n), eltype(cα), eltype(c))
    return atomcoordinates(bp, dihedrals, SVector{3,T}(n), SVector{3,T}(cα), SVector{3,T}(c))
end
atomcoordinates(bp::BondParametrization, dihedrals::Vector, n::Atom, cα::Atom, c::Atom) = atomcoordinates(bp, dihedrals, n.coords, cα.coords, c.coords)
function atomcoordinates(bp::BondParametrization, dihedrals::Vector, chain::Chain)
    nter = first(chain)::Residue
    return atomcoordinates(bp, dihedrals, nter["N"]::Atom, nter["CA"]::Atom, nter["C"]::Atom)
end

"""
    out = buildchain(reference::Chain, bp::BondParametrization, X::AbstractVector{<:SVector{3}})

Copy `reference` and replace its atom coordinates with `X`, matched through
`bp.atoms`.
"""
function buildchain(reference::Chain, bp::BondParametrization, X::AbstractVector{<:SVector{3}})
    out = copy(reference)
    coordidx = Dict{AtomData, Int}()
    for (i, adata) in enumerate(bp.atoms)
        coordidx[adata] = i
    end
    for a in collectatoms(out)
        i = coordidx[AtomData(a)]
        a.coords .= X[i]
    end
    return out
end
