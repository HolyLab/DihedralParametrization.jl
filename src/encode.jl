# Encoding is optimized for repeated calls to `atomcoordinates`.

"""
    AtomData

Identifies one atom of a `BondParametrization` by its source residue and
atom name.

# Fields
- `ridx::Int`: the residue number from the source structure (as returned by
  `BioStructures.resnumber`), not a position in `BondParametrization.residues`.
- `aname::Symbol`: the atom's name within its residue (e.g. `:CA`, `:HB2`).
"""
struct AtomData
    ridx::Int
    aname::Symbol
end
function AtomData(a::Atom)
    aname = atomname(a)
    return AtomData(parse(Int, resid(residue(a); full=false)), Symbol(aname))
end

# A build step that places one atom `d` by extending the chain `a-b-c-d` with
# the SN-NeRF formula (see `snnerf`): `predecessors` are the indices into the
# coordinate vector of the already-placed atoms `a`, `b`, `c`; `ℓcd` and
# `θbcd` are the fixed C–D bond length and B–C–D bond angle; `ϕ` is the
# dihedral angle to use when `rotatable` is `false` (a rotatable step instead
# consumes the next entry of `dihedrals`).
struct Extend{T}
    predecessors::Tuple{Int,Int,Int}              # indices in X of a, b, c
    ℓcd::T
    θbcd::T
    rotatable::Bool
    ϕ::T  # fixed (or original) dihedral angle, only used if not rotatable
end

# A build step that places one or more atoms in a fixed (non-rotatable)
# tetrahedral-style geometry from already-placed atoms `a`, `b`, `c` (indexed
# by `predecessors`), via `add_to_middle!`. `βs` holds each new atom's
# coefficients in the local `a-b-c` frame.
struct Branch{T}
    predecessors::Tuple{Int,Int,Int}              # indices in X of a, b, c
    βs::Vector{SVector{3,T}}                      # coefficients for placement from a, b, c
end

# One residue's sidechain build sequence: the ordered list of `Extend`/`Branch`
# steps that places its atoms once the backbone and preceding residues are
# already placed.
struct ResidueData{T}
    steps::Vector{Union{Extend{T}, Branch{T}}}    # list of build steps for this residue
end


"""
    BondParametrization{T<:Real}

Stores the fixed geometry and build sequence needed to reconstruct a protein
chain from its rotatable dihedral angles. See `bondparametrization`,
`atomcoordinates`, and `buildchain`. `ndihedrals(bp)` gives the required
number of angles.

# Fields
- `atoms::Vector{AtomData}`: one entry per atom, in the order
  `atomcoordinates` returns coordinates and `buildchain` consumes them; entry
  `i` names the atom at row `i` of that coordinate vector.
- `bblengths::Vector{T}`: backbone bond lengths, length `3*nres` where
  `nres = length(residues)`. Entries `3i-2:3i` for residue `i` are the N–Cα,
  Cα–C, and C–N(next) lengths; the last residue's "next" atom is OXT.
- `bbangles::Vector{T}`: backbone bond angles, length `3*nres-1`. Entries
  `3i-2:3i` for residue `i` are the N–Cα–C, Cα–C–N(next), and
  C–N(next)–Cα(next) angles; the last residue contributes only the first two
  of these (there is no following residue).
- `omegas::Vector{T}`: the ω dihedral of each peptide bond, length
  `nres - 1`. These are fixed and do not appear in `dihedrals`; entry `i` is
  the ω between residue `i` and residue `i+1`.
- `phirotatable::Vector{Bool}`: whether residue `i`'s own φ is a free
  dihedral (`true`) or fixed by a ring closure such as proline's (`false`),
  length `nres`. Residue 1 has no φ; entry 1 is unused.
- `phi::Vector{T}`: the fixed φ value for residue `i`, length `nres`; only
  meaningful where `!phirotatable[i]`.
- `residues::Vector{ResidueData{T}}`: each residue's sidechain build
  sequence — the bond length, bond angle, and rotatability of every
  sidechain atom placed relative to its already-placed predecessors.
"""
struct BondParametrization{T<:Real}
    atoms::Vector{AtomData}    # list of all atoms in the chain
    bblengths::Vector{T}          # backbone bond lengths, of length 3*nres (we use OXT for the last C)
    bbangles::Vector{T}           # backbone bond angles, of length 3*nres-1 (")
    omegas::Vector{T}             # omega dihedrals (fixed and not represented in the dihedrals vector), of length nres-1
    phirotatable::Vector{Bool}    # whether residue i's own phi is a free dihedral, of length nres (residue 1 has no phi; entry unused)
    phi::Vector{T}                # fixed phi value for residue i, of length nres (only meaningful where !phirotatable[i])
    residues::Vector{ResidueData{T}}   # list of residues
end

Base.show(io::IO, bp::BondParametrization) = print(io, "BondParametrization with $(length(bp.atoms)) atoms and $(length(bp.residues)) residues")

"""
    n = ndihedrals(bp::BondParametrization)

Return the number of rotatable dihedrals in `bp`.
"""
function ndihedrals(bp::BondParametrization)
    nres = length(bp.residues)
    nphi_free = count(view(bp.phirotatable, 2:nres))
    nsidechain_free = sum((count(step -> step isa Extend && step.rotatable, r.steps) for r in bp.residues); init=0)
    return (nres - 1) + nphi_free + 1 + nsidechain_free
end

# Error for a residue absent from the build tables.
unrecognized_residue_message(i::Int, rname::AbstractString) =
    "residue $i ($rname): unrecognized residue name — " *
    "histidine must be disambiguated as HID/HIE/HIP, and amino-/carboxyl-terminal " *
    "residues need an \"N\"/\"C\"-prefixed name; BioStructures.specializeresnames! " *
    "assigns these names (see the documentation's \"Pre-requisites\" section)"

"""
    resatom, ridx = resolvebuildref(ref::AbstractString, i::Int, ress, resatomidxs)

Resolve a build-step reference atom name `ref` for residue `ress[i]`.
Standard names resolve within-residue. A name prefixed with `-` or `+`
(e.g. `"-C"`, `"+N"`, the CHARMM convention) names an atom in the previous
or next residue rather than in `ress[i]`.

Returns the `Atom` together with its index into the coordinate array that
`atomcoordinates` builds.
"""
function resolvebuildref(ref::AbstractString, i::Int, ress, resatomidxs)
    c1 = ref[1]
    if c1 == '-' || c1 == '+'
        aname = ref[2:end]
        j = c1 == '-' ? i - 1 : i + 1
        ok = 1 <= j <= length(ress) && (c1 == '-' ? sequentialresidues(ress[j], ress[i]) : sequentialresidues(ress[i], ress[j]))
        ok || error("residue $i ($(resname(ress[i]))): cannot resolve build-step reference \"$ref\" — " *
                    (1 <= j <= length(ress) ? "residue $j is not sequential with residue $i (chain break)" :
                     "residue $i has no " * (c1 == '-' ? "preceding" : "following") * " residue"))
        return ress[j][aname], resatomidxs[j][aname]
    else
        atom = try
            ress[i][ref]
        catch err
            err isa KeyError || rethrow()
            error("residue $i ($(resname(ress[i]))): cannot find atom \"$ref\" required by residue_build_sequence")
        end
        return atom, resatomidxs[i][ref]
    end
end

"""
    bp, dihedrals = bondparametrization(chain)
    bp, dihedrals = bondparametrization(T, chain)

Represent a protein `chain` by fixed bond parameters `bp` and its rotatable
`dihedrals`. These values determine its coordinates up to a rigid transform.

`dihedrals` is a vector representing the rotatable dihedral angles (in
radians) in the chain, in the order documented in [Conventions](@ref);
`ndihedrals` gives its length. `bp` is a `BondParametrization` object
containing all other necessary information.

The first form stores `bp`'s bond lengths, bond angles, and fixed dihedrals
(and returns `dihedrals`) as `Float64`; pass a type explicitly, e.g.
`bondparametrization(Float32, chain)`, to use a different real type instead.
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
    phirotatable = Vector{Bool}(undef, nres)
    phi = Vector{T}(undef, nres)
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
            # A ring through N-CA fixes ϕ (as in proline).
            rot = try
                residue_phi_rotatable[resname(ress[i+1])]
            catch err
                err isa KeyError || rethrow()
                error(unrecognized_residue_message(i+1, resname(ress[i+1])))
            end
            phirotatable[i+1] = rot
            phi[i+1] = ϕs[i+1]
            rot && push!(dihedrals, ϕs[i+1])
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
        # CYX is disulfide-bonded cysteine and therefore has no thiol hydrogen.
        rname ∈ ("CYX", "NCYX", "CCYX") && haskey(res.atoms, "HG") &&
            error("residue $i ($rname): CYX (disulfide-bonded cysteine) must not have a thiol HG")
        seq = try
            residue_build_sequence[rname]
        catch err
            err isa KeyError || rethrow()
            error(unrecognized_residue_message(i, rname))
        end
        atomidx = resatomidxs[i]
        steps = Vector{Union{Extend{T}, Branch{T}}}()
        for step in seq
            if isa(step, Tuple{String,String,String,String,Bool})
                # Extend
                a, b, c, d, rotatable = step
                i == length(ress) && d == "OXT" && continue # already added
                atomA, idxA = resolvebuildref(a, i, ress, resatomidxs)
                atomB, idxB = resolvebuildref(b, i, ress, resatomidxs)
                atomC, idxC = resolvebuildref(c, i, ress, resatomidxs)
                predecessors = (idxA, idxB, idxC)
                atomD = try
                    res[d]
                catch err
                    err isa KeyError || rethrow()
                    error("residue $i ($rname): cannot find atom \"$d\" required by residue_build_sequence")
                end
                ℓcd = norm(atomC.coords - atomD.coords)
                θbcd = bondangle(atomB, atomC, atomD)
                ϕ = dihedralangle(atomA, atomB, atomC, atomD)
                push!(steps, Extend{T}(predecessors, ℓcd, θbcd, rotatable, ϕ))
                if rotatable
                    push!(dihedrals, ϕ)
                end
                addatom(atomD)
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
    return BondParametrization{T}(atoms, bblengths, bbangles, omegas, phirotatable, phi, residues), dihedrals
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
