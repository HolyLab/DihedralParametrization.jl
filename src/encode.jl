# Note: encoding happens once, whereas decoding happens repeatedly during optimization.
# Therefore, it is acceptable for encoding to be slower and more complex if it speeds up decoding.
# Thus, the implementation here is optimized to make `atomcoordinates` as fast as possible.

struct AtomData
    ridx::Int
    aname::Symbol
end
function AtomData(a::Atom)
    aname = atomname(a)
    return AtomData(parse(Int, resid(residue(a); full=false)), Symbol(aname))
end

struct Extend{T}
    predecessors::Tuple{Int,Int,Int}              # indices in X of a, b, c
    ℓcd::T
    θbcd::T
    rotatable::Bool
    ϕ::T  # fixed (or original) dihedral angle, only used if not rotatable
end

struct Branch{T}
    predecessors::Tuple{Int,Int,Int}              # indices in X of a, b, c
    βs::Vector{SVector{3,T}}                      # coefficients for placement from a, b, c
end

struct ResidueData{T}
    steps::Vector{Union{Extend{T}, Branch{T}}}    # list of build steps for this residue
end


struct BondParametrization{T<:Real}
    atoms::Vector{AtomData}    # list of all atoms in the chain
    bblengths::Vector{T}          # backbone bond lengths, of length 3*nres (we use OXT for the last C)
    bbangles::Vector{T}           # backbone bond angles, of length 3*nres-1 (")
    omegas::Vector{T}             # omega dihedrals (fixed and not represented in the dihedrals vector), of length nres (final one is for placement of OXT)
    residues::Vector{ResidueData{T}}   # list of residues
end

Base.show(io::IO, bp::BondParametrization) = print(io, "BondParametrization with $(length(bp.atoms)) atoms and $(length(bp.residues)) residues")

"""
    bp, dihedrals = bondparametrization(chain)

Parametrize a protein `chain` via bond parameters. Together, `bp` and
`dihedrals` are sufficient to reconstruct the 3D coordinates of the chain modulo
a rigid body transformation.

`dihedrals` is a vector representing the rotatable dihedral angles in the chain.
`bp` is a `BondParametrization` object containing all other necessary information.
"""
function bondparametrization(::Type{T}, chain::Chain) where {T<:Real}
    # Preallocate arrays
    ress = collectresidues(chain)
    nres = length(ress)
    natoms = sum(length(res) for res in ress)
    atoms = Vector{AtomData}(undef, natoms)
    bblengths = Vector{T}(undef, 3nres)
    bbangles = Vector{T}(undef, 3nres-1)
    omegas = deleteat!(omegaangles(ress), 1)   # we'll add one more omega at the end for placement of OXT
    residues = Vector{ResidueData{T}}(undef, nres)
    dihedrals = T[]

    atomidx = Dict{String,Int}()
    resatomidxs = Vector{typeof(atomidx)}(undef, nres)
    aidx = 0
    function addatom(a::Atom)
        aidx += 1
        atoms[aidx] = AtomData(a)
        atomidx[atomname(a)] = aidx
    end

    # Backbone atoms
    ϕs, ψs = phiangles(ress), psiangles(ress)
    for (i, res) in enumerate(ress)
        empty!(atomidx)
        n, cα, c = res["N"], res["CA"], res["C"]
        addatom(n); addatom(cα); addatom(c)
        next = if i < nres
            ress[i+1]["N"]
        else
            a = res["OXT"]
            addatom(a)
            a
        end
        bblengths[3i - 2] = norm(n.coords - cα.coords)      # N - Cα
        bblengths[3i - 1] = norm(cα.coords - c.coords)      # Cα - C
        bblengths[3i]     = norm(c.coords - next.coords)    # C - N(next)
        bbangles[3i - 2] = bondangle(n, cα, c)     # N - Cα - C
        bbangles[3i - 1] = bondangle(cα, c, next)  # Cα - C - N(next)
        if i < nres
            push!(dihedrals, ψs[i])
            push!(dihedrals, ϕs[i+1])
            nextnext = ress[i+1]["CA"]
            bbangles[3i] = bondangle(c, next, nextnext)  # C - N(next) - Cα(next)
        else
            push!(dihedrals, dihedralangle(n, cα, c, next))  # use OXT for setting the final ψ
        end
        resatomidxs[i] = copy(atomidx)
    end
    # Sidechain atoms
    for (i, res) in enumerate(ress)
        rname = resname(res)
        seq = residue_build_sequence[rname]
        atomidx = resatomidxs[i]
        steps = Vector{Union{Extend{T}, Branch{T}}}()
        for step in seq
            if isa(step, Tuple{String,String,String,String,Bool})
                # Extend
                a, b, c, d, rotatable = step
                i == length(ress) && d == "OXT" && continue # already added
                if i == 1 && d == "H" && !haskey(res.atoms, "H")  # FIXME: this is a hack
                    d = "H2"
                end
                predecessors = (atomidx[a], atomidx[b], atomidx[c])
                ℓcd = norm(res[c].coords - res[d].coords)
                θbcd = bondangle(res[b], res[c], res[d])
                ϕ = dihedralangle(res[a], res[b], res[c], res[d])
                push!(steps, Extend{T}(predecessors, ℓcd, θbcd, rotatable, ϕ))
                if rotatable
                    push!(dihedrals, ϕ)
                end
                addatom(res[d])
            elseif isa(step, Tuple{String,String,String,Vector{String}})
                # Branch
                a, b, c, ats = step
                i == length(ress) && length(ats) == 1 && only(ats) == "OXT" && continue # already added
                predecessors = (atomidx[a], atomidx[b], atomidx[c])
                βs = betas(a, b, c, ats, res)
                push!(steps, Branch{T}(predecessors, βs))
                for at in ats
                    addatom(res[at])
                end
            else
                error("Invalid step: $step")
            end
        end
        residues[i] = ResidueData{T}(steps)
    end
    @assert aidx == natoms   # amino terminus may have two extra H
    return BondParametrization{T}(atoms, bblengths, bbangles, omegas, residues), dihedrals
end
bondparametrization(chain::Chain) = bondparametrization(Float64, chain)

betas(a::AbstractString, b::AbstractString, c::AbstractString, ats, res) = betas(SVector{3}(res[a].coords), SVector{3}(res[b].coords), SVector{3}(res[c].coords), [SVector{3}(res[d].coords) for d in ats])
function betas(a::AbstractVector{T}, b::AbstractVector{T}, c::AbstractVector{T}, ds) where T<:Real
    # The mirror image of add_to_middle
    ab = b - a
    ab = ab / norm(ab)
    cb = b - c
    cb = cb / norm(cb)
    cb = cb / norm(cb)
    n = cross(ab, cb)
    n = n / norm(n)
    Minv = inv([ab cb n])
    return map(ds) do d
        db = d - b
        Minv * db
    end
end
