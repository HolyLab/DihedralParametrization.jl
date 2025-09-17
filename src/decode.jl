"""
    d = snnerf(a, b, c, ℓcd, θbcd, ϕbc)

Given three points `a`, `b`, and `c` (as 3D vectors), the bond length `ℓcd`,
the bond angle `θbcd` (in radians), and the dihedral angle `ϕbc` (in radians),
compute the coordinates of point `d` using the SN-NeRF algorithm.

## Reference

> Parsons, Jerod, et al. "Practical conversion from torsion space to Cartesian
> space for in silico protein synthesis." Journal of computational chemistry 26.10 (2005): 1063-1068.
"""
function snnerf(a::AbstractVector{T}, b::AbstractVector{T}, c::AbstractVector{T},
                ℓcd::T, θbcd::T, ϕbc::T) where T<:Real
    return snnerf(a, b, c, ℓcd, sincos(θbcd), sincos(ϕbc))
end

function snnerf(a::AbstractVector{T}, b::AbstractVector{T}, c::AbstractVector{T},
                ℓcd::T, θbcd::T, (sϕ, cϕ)::Tuple{T,T}) where T<:Real
    return snnerf(a, b, c, ℓcd, sincos(θbcd), (sϕ, cϕ))
end

function snnerf(a::AbstractVector{T}, b::AbstractVector{T}, c::AbstractVector{T},
                ℓcd::T, (sθ, cθ)::Tuple{T,T}, (sϕ, cϕ)::Tuple{T,T}) where T<:Real
    # Calculate unit vectors
    bc = c - b
    bc = bc / norm(bc)    # technically, SN-NeRF wants us to pass `norm(bc)` in as an argument
    ab = b - a
    n = cross(ab, bc)
    n = n / norm(n)
    M = [bc cross(n, bc) n]
    ℓsθ, ℓcθ = ℓcd .* (sθ, cθ)
    d2 = SVector(-ℓcθ, ℓsθ * cϕ, ℓsθ * sϕ)
    return c + M * d2
end

function add_to_middle(a::AbstractVector{T}, b::AbstractVector{T}, c::AbstractVector{T}, βs::AbstractVector{T}...) where T<:Real
    # used, e.g., to place the Cβ atom in a tetrahedral geometry
    ab = b - a
    ab = ab / norm(ab)
    cb = b - c
    cb = cb / norm(cb)
    n = cross(ab, cb)
    n = n / norm(n)
    M = [ab cb n]
    return Ref(b) .+ Ref(M) .* βs
end

"""
    X = atomcoordinates(bp::BondParametrization{T}, dihedrals::Vector{T}, n::AbstractVector{T}, cα::AbstractVector{T}, c::AbstractVector{T}) where T<:Real

Given a `BondParametrization` object `bp`, a vector of dihedral angles `dihedrals`,
and the coordinates of the first three backbone atoms in the chain (`n`, `cα`, and `c`),
compute the 3D coordinates of all atoms in the chain.
"""
function atomcoordinates(bp::BondParametrization, dihedrals::Vector{T}, n::SVector{3,T}, cα::SVector{3,T}, c::SVector{3,T}) where {T<:Real}
    # Check that the inputs are consistent with `bp`
    norm(cα - n) ≈ bp.bblengths[1] || error("Provided N and Cα do not match bond length in bp")
    norm(c - cα) ≈ bp.bblengths[2] || error("Provided Cα and C do not match bond length in bp")
    bondangle(n - cα, c - cα) ≈ bp.bbangles[1] || error("Provided N, Cα, and C do not match bond angle in bp")

    # Initialize the coordinates array
    X = sizehint!([n, cα, c], length(bp.atoms))
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
        φ = dihedrals[idx+=1]       # rotatable
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
                d = add_to_middle(a, b, c, step.βs...)
                push!(X, d...)
            end
        end
    end
    @assert length(X) == length(bp.atoms)
    return X
end
atomcoordinates(bp::BondParametrization, dihedrals::Vector, n::Atom, cα::Atom, c::Atom) = atomcoordinates(bp, dihedrals, SVector{3}(n.coords), SVector{3}(cα.coords), SVector{3}(c.coords))
function atomcoordinates(bp::BondParametrization, dihedrals::Vector, chain::Chain)
    nter = first(chain)
    return atomcoordinates(bp, dihedrals, nter["N"], nter["CA"], nter["C"])
end

"""
    out = buildchain(reference::Chain, bp::BondParametrization, X::AbstractVector{<:SVector{3}})

Given a reference `Chain` object, a `BondParametrization` object `bp`, and a
vector of atom 3D coordinates `X`, construct a new `Chain` object with the same
sequence and atoms as `reference` but with the coordinates from `X`.
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
