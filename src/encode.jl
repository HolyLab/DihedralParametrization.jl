# Encoding is optimized for repeated calls to `atomcoordinates`.

"""
    AtomData

Identifies one atom of a `BondParametrization` by its source residue and
atom name.

# Fields
- `ridx::Int`: the residue number from the source structure (as returned by
  `BioStructures.resnumber`), not a position in `BondParametrization.atoms`.
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
# coordinate vector of the already-placed atoms `a`, `b`, `c`, and `aidx` is
# the index `d` itself receives. `ℓcd` and `θbcd` are the fixed C–D bond
# length and B–C–D bond angle; `ϕ` is the dihedral angle to use when
# `rotatable` is `false` (a rotatable step instead consumes the next entry of
# `dihedrals`).
struct Extend{T}
    predecessors::Tuple{Int,Int,Int}              # indices in X of a, b, c
    aidx::Int                                     # index in X of the atom this step places
    ℓcd::T
    θbcd::T
    rotatable::Bool
    ϕ::T  # fixed (or original) dihedral angle, only used if not rotatable
end

# A build step that places one or more atoms in a fixed (non-rotatable)
# tetrahedral-style geometry from already-placed atoms `a`, `b`, `c` (indexed
# by `predecessors`), via `add_to_middle!`. `βs` holds each new atom's
# coefficients in the local `a-b-c` frame; the atoms occupy indices
# `aidx:aidx+length(βs)-1`.
struct Branch{T}
    predecessors::Tuple{Int,Int,Int}              # indices in X of a, b, c
    aidx::Int                                     # index in X of the first atom this step places
    βs::Vector{SVector{3,T}}                      # coefficients for placement from a, b, c
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
- `nres::Int`: the number of residues in the chain.
- `ℓnca::T`, `ℓcac::T`, `θncac::T`: the N–Cα bond length, Cα–C bond length,
  and N–Cα–C bond angle of residue 1. These fix the reference frame;
  `atomcoordinates` validates the N, Cα, and C coordinates it is given
  against them.
- `steps::Vector{Union{Extend{T},Branch{T}}}`: the build sequence. Atoms
  `1:3` are residue 1's N, Cα, and C; every later atom is placed by exactly
  one step, and the steps appear in the order their atoms appear in `atoms`.
  Each step names the already-placed atoms it builds from, the bond length
  and bond angle (or, for a `Branch`, the local-frame coefficients) that fix
  its geometry, and whether its dihedral is rotatable.
- `ndihedrals::Int`: the number of rotatable steps, and hence the required
  length of `dihedrals`.
"""
struct BondParametrization{T<:Real}
    atoms::Vector{AtomData}
    nres::Int
    ℓnca::T                       # residue 1: N–Cα bond length
    ℓcac::T                       # residue 1: Cα–C bond length
    θncac::T                      # residue 1: N–Cα–C bond angle
    steps::Vector{Union{Extend{T},Branch{T}}}
    ndihedrals::Int
end

Base.show(io::IO, bp::BondParametrization) = print(io, "BondParametrization with $(length(bp.atoms)) atoms and $(bp.nres) residues")

"""
    n = ndihedrals(bp::BondParametrization)

Return the number of rotatable dihedrals in `bp`.
"""
ndihedrals(bp::BondParametrization) = bp.ndihedrals

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
    ress = collectresidues(chain)
    nres = length(ress)
    nres >= 1 || error("chain has no residues")
    natoms = sum(length(res) for res in ress)
    atoms = Vector{AtomData}(undef, natoms)
    steps = Vector{Union{Extend{T},Branch{T}}}()
    dihedrals = T[]

    atomidx = Dict{String,Int}()
    resatomidxs = Vector{typeof(atomidx)}(undef, nres)
    aidx = 0
    function addatom(a::Atom)
        aidx += 1
        atoms[aidx] = AtomData(a)
        atomidx[atomname(a)] = aidx
        return aidx
    end

    # Backbone atoms. Residue 1's N, Cα, and C are the reference frame; every
    # later backbone atom gets a step.
    ϕs, ψs, ωs = phiangles(ress), psiangles(ress), omegaangles(ress)
    ℓnca = ℓcac = θncac = zero(T)   # residue 1's frame, measured on the first iteration
    previdx = (0, 0, 0)
    for (i, res) in enumerate(ress)
        empty!(atomidx)
        n, cα, c = res["N"], res["CA"], res["C"]
        idxn, idxcα, idxc = addatom(n), addatom(cα), addatom(c)
        if i == 1
            ℓnca = norm(n.coords - cα.coords)
            ℓcac = norm(cα.coords - c.coords)
            θncac = bondangle(n, cα, c)
        else
            prevcα, prevc = ress[i-1]["CA"], ress[i-1]["C"]
            # N_i rotates about the Cα_{i-1}–C_{i-1} bond: ψ_{i-1}.
            push!(steps, Extend{T}(previdx, idxn, norm(prevc.coords - n.coords),
                                   bondangle(prevcα, prevc, n), true, ψs[i-1]))
            push!(dihedrals, ψs[i-1])
            # Cα_i is fixed by the peptide bond's ω.
            push!(steps, Extend{T}((previdx[2], previdx[3], idxn), idxcα, norm(n.coords - cα.coords),
                                   bondangle(prevc, n, cα), false, ωs[i]))
            # C_i rotates about the N_i–Cα_i bond: φ_i. A ring through N–Cα
            # fixes φ (as in proline).
            rot = try
                residue_phi_rotatable[resname(res)]
            catch err
                err isa KeyError || rethrow()
                error(unrecognized_residue_message(i, resname(res)))
            end
            push!(steps, Extend{T}((previdx[3], idxn, idxcα), idxc, norm(cα.coords - c.coords),
                                   bondangle(n, cα, c), rot, ϕs[i]))
            rot && push!(dihedrals, ϕs[i])
        end
        if i == nres
            # OXT closes the carboxyl terminus; its dihedral is a final ψ.
            oxt = res["OXT"]
            idxoxt = addatom(oxt)
            ψ′ = dihedralangle(n, cα, c, oxt)
            push!(steps, Extend{T}((idxn, idxcα, idxc), idxoxt, norm(c.coords - oxt.coords),
                                   bondangle(cα, c, oxt), true, ψ′))
            push!(dihedrals, ψ′)
        end
        previdx = (idxn, idxcα, idxc)
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
        for step in seq
            if isa(step, Tuple{String,String,String,String,Bool})
                # Extend
                a, b, c, d, rotatable = step
                i == nres && d == "OXT" && continue # already added
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
                push!(steps, Extend{T}(predecessors, addatom(atomD), ℓcd, θbcd, rotatable, ϕ))
                if rotatable
                    push!(dihedrals, ϕ)
                end
            elseif isa(step, Tuple{String,String,String,Vector{String}})
                # Branch
                a, b, c, ats = step
                i == nres && length(ats) == 1 && only(ats) == "OXT" && continue # already added
                predecessors = (atomidx[a], atomidx[b], atomidx[c])
                βs = betas(a, b, c, ats, res)
                push!(steps, Branch{T}(predecessors, aidx + 1, βs))
                for at in ats
                    addatom(res[at])
                end
            else
                error("Invalid step: $step")
            end
        end
    end
    if aidx != natoms
        placed = Set(view(atoms, 1:aidx))
        detail = ""
        for (i, res) in enumerate(ress)
            unplaced = sort!([name for name in keys(res.atoms) if AtomData(res[name]) ∉ placed])
            if !isempty(unplaced)
                detail = "; residue $i ($(resname(res))) has atom(s) " * join(unplaced, ", ") *
                         " that residue_build_sequence never places"
                break
            end
        end
        error("the build sequence placed $aidx atoms but the chain has $natoms$detail")
    end
    ndih = count(step -> step isa Extend && step.rotatable, steps)
    length(dihedrals) == ndih ||
        error("recorded $(length(dihedrals)) dihedrals for $ndih rotatable steps")
    return BondParametrization{T}(atoms, nres, ℓnca, ℓcac, θncac, steps, ndih), dihedrals
end
bondparametrization(chain::Chain) = bondparametrization(Float64, chain)

betas(a::AbstractString, b::AbstractString, c::AbstractString, ats, res) = betas(SVector{3}(res[a].coords), SVector{3}(res[b].coords), SVector{3}(res[c].coords), [SVector{3}(res[d].coords) for d in ats])
function betas(a::AbstractVector{T}, b::AbstractVector{T}, c::AbstractVector{T}, ds) where T<:Real
    # The mirror image of add_to_middle
    ab = b - a
    ab = ab / norm(ab)
    cb = b - c
    cb = cb / norm(cb)
    n = cross(ab, cb)
    n = n / norm(n)
    Minv = inv([ab cb n])
    return map(ds) do d
        db = d - b
        Minv * db
    end
end
