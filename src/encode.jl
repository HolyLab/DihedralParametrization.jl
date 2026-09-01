"""
    AtomKey

Identifies one atom of a `BondParametrization` by residue number and atom
name.

# Fields
- `resnum::Int`: the residue number from the source structure (as returned by
  `BioStructures.resnumber`), not a position in `BondParametrization.atoms`.
- `aname::Symbol`: the atom's name within its residue (e.g. `:CA`, `:HB2`).
"""
struct AtomKey
    resnum::Int
    aname::Symbol
end
function AtomKey(a::Atom)
    res = residue(a)
    # An insertion code would make the residue number ambiguous as a key.
    inscode(res) == ' ' ||
        throw(ArgumentError("$(resdesc(res)): insertion codes are not supported"))
    return AtomKey(resnumber(res), Symbol(atomname(a)))
end

# Name a residue in an error message by its number in the source structure
# (with its insertion code, if it has one) and its name, never by its
# position in the chain.
resdesc(res::Residue) = "residue $(resid(res; full=false)) ($(resname(res)))"

# Fetch a backbone atom (N, CA, C, or the carboxyl terminus's OXT) from `res`.
function backboneatom(res::Residue, name::AbstractString)
    haskey(res.atoms, name) ||
        throw(ArgumentError("$(resdesc(res)): missing backbone atom \"$name\""))
    return res[name]
end

# Fetch an atom that `residue_build_sequence` names for `res`.
function templateatom(res::Residue, name::AbstractString)
    haskey(res.atoms, name) ||
        throw(ArgumentError("$(resdesc(res)): cannot find atom \"$name\" required by residue_build_sequence"))
    return res[name]
end

# Compare and hash these mutable-field values by contents.
_fieldsmatch(eq, x, y) = all(f -> eq(getfield(x, f), getfield(y, f)), fieldnames(typeof(x)))
_hashfields(x, h::UInt) = foldr((f, h) -> hash(getfield(x, f), h), fieldnames(typeof(x)); init=h)

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

Base.:(==)(x::Extend, y::Extend) = _fieldsmatch(==, x, y)
Base.isequal(x::Extend, y::Extend) = _fieldsmatch(isequal, x, y)
Base.hash(x::Extend, h::UInt) = _hashfields(x, hash(:Extend, h))

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

Base.:(==)(x::Branch, y::Branch) = _fieldsmatch(==, x, y)
Base.isequal(x::Branch, y::Branch) = _fieldsmatch(isequal, x, y)
Base.hash(x::Branch, h::UInt) = _hashfields(x, hash(:Branch, h))

"""
    BondParametrization{T<:Real}

Fixed geometry and build sequence for reconstructing a protein chain from
its rotatable dihedral angles. See `bondparametrization`, `atomcoordinates`,
and `buildchain`.

# Fields
- `atoms::Vector{AtomKey}`: one entry per atom, in the order
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

Use [`natoms`](@ref) and [`ndihedrals`](@ref) instead of accessing fields.
"""
struct BondParametrization{T<:Real}
    atoms::Vector{AtomKey}
    nres::Int
    ℓnca::T                       # residue 1: N–Cα bond length
    ℓcac::T                       # residue 1: Cα–C bond length
    θncac::T                      # residue 1: N–Cα–C bond angle
    steps::Vector{Union{Extend{T},Branch{T}}}
    ndihedrals::Int

    function BondParametrization{T}(atoms, nres, ℓnca, ℓcac, θncac, steps, ndihedrals) where {T<:Real}
        nrot = count(step -> step isa Extend && step.rotatable, steps)
        nrot == ndihedrals ||
            throw(ArgumentError("steps has $nrot rotatable steps but ndihedrals = $ndihedrals"))
        return new{T}(atoms, nres, ℓnca, ℓcac, θncac, steps, ndihedrals)
    end
end

Base.show(io::IO, bp::BondParametrization) = print(io, "BondParametrization with $(length(bp.atoms)) atoms and $(bp.nres) residues")

Base.:(==)(x::BondParametrization, y::BondParametrization) = _fieldsmatch(==, x, y)
Base.isequal(x::BondParametrization, y::BondParametrization) = _fieldsmatch(isequal, x, y)
Base.hash(x::BondParametrization, h::UInt) = _hashfields(x, hash(:BondParametrization, h))

Base.eltype(::Type{BondParametrization{T}}) where {T} = T

"""
    n = ndihedrals(bp::BondParametrization)
    n = ndihedrals(plan::JacobianPlan)

Return the number of rotatable dihedrals.
"""
ndihedrals(bp::BondParametrization) = bp.ndihedrals

"""
    n = natoms(bp::BondParametrization)
    n = natoms(plan::JacobianPlan)

Return the number of atoms.
"""
natoms(bp::BondParametrization) = length(bp.atoms)

# Error for a residue absent from the build tables.
unrecognized_residue_message(res::Residue) =
    "$(resdesc(res)): unrecognized residue name — " *
    "histidine must be disambiguated as HID/HIE/HIP, and amino-/carboxyl-terminal " *
    "residues need an \"N\"/\"C\"-prefixed name; BioStructures.specializeresnames! " *
    "assigns these names (see the documentation's \"Pre-requisites\" section)"

"""
    resatom, aidx = resolvebuildref(ref::AbstractString, i::Int, ress, resatomidxs)

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
        ok || throw(ArgumentError("$(resdesc(ress[i])): cannot resolve build-step reference \"$ref\" — " *
                                  (1 <= j <= length(ress) ?
                                   "residue $(resnumber(ress[j])) is not sequential with residue $(resnumber(ress[i])) (chain break)" :
                                   "it has no " * (c1 == '-' ? "preceding" : "following") * " residue")))
        return templateatom(ress[j], aname), resatomidxs[j][aname]
    else
        return templateatom(ress[i], ref), resatomidxs[i][ref]
    end
end

"""
    dihedrals = dihedralangles(bp::BondParametrization, X::AbstractVector{<:SVector{3}})
    dihedrals = dihedralangles(bp::BondParametrization, chain::Chain)

Measure a configuration's rotatable dihedral angles in radians. This is the
inverse of `atomcoordinates`; `dihedrallabels` names the entries. Angles lie
in `-π` to `π`.

`X` holds one coordinate per entry of `bp.atoms`, as `atomcoordinates`
returns; a `DimensionMismatch` is thrown if its length differs. The `Chain`
method takes the coordinates from `chain`, matched through `bp.atoms`, and
throws an `ArgumentError` if the chain's atom count differs from `bp`'s or if
an atom named in `bp.atoms` is absent.

The result promotes the element types of `bp` and the coordinates. See
[`dihedralangles!`](@ref) for the in-place form.
"""
function dihedralangles(bp::BondParametrization, X::AbstractVector{<:SVector{3}})
    R = promote_type(eltype(bp), eltype(eltype(X)))
    return dihedralangles!(Vector{R}(undef, ndihedrals(bp)), bp, X)
end

function dihedralangles(bp::BondParametrization, chain::Chain)
    return dihedralangles(bp, chaincoordinates(bp, chain))
end

"""
    dihedralangles!(dihedrals, bp::BondParametrization, X::AbstractVector{<:SVector{3}})

Write [`dihedralangles`](@ref) into `dihedrals` and return it. The output must
have `ndihedrals(bp)` entries.
"""
function dihedralangles!(dihedrals::AbstractVector, bp::BondParametrization, X::AbstractVector{<:SVector{3}})
    Base.require_one_based_indexing(X)
    length(X) == length(bp.atoms) ||
        throw(DimensionMismatch("length(X) = $(length(X)) does not match bp's $(length(bp.atoms)) atoms"))
    nd = ndihedrals(bp)
    length(dihedrals) == nd ||
        throw(DimensionMismatch("length(dihedrals) = $(length(dihedrals)) does not match bp's $nd rotatable dihedrals"))
    k = firstindex(dihedrals) - 1
    for step in bp.steps
        (step isa Extend && step.rotatable) || continue
        a, b, c = step.predecessors
        dihedrals[k+=1] = dihedralangle(X[b] - X[a], X[c] - X[b], X[step.aidx] - X[c])
    end
    return dihedrals
end

# Collect `chain`'s coordinates in the order of `bp.atoms`.
function chaincoordinates(bp::BondParametrization, chain::Chain)
    chainatoms = collectatoms(chain)
    length(chainatoms) == length(bp.atoms) ||
        throw(ArgumentError("chain has $(length(chainatoms)) atoms but the parametrization describes $(length(bp.atoms))"))
    byatomkey = Dict{AtomKey,eltype(chainatoms)}()
    for a in chainatoms
        byatomkey[AtomKey(a)] = a
    end
    return map(bp.atoms) do akey
        a = get(byatomkey, akey, nothing)
        a === nothing &&
            throw(ArgumentError("chain has no atom $(akey.aname) in residue $(akey.resnum), which the parametrization requires"))
        SVector{3}(a.coords)
    end
end

"""
    DihedralLabel

Names one entry of a rotatable-dihedral vector by the atoms that define it.
See `dihedrallabels`.

# Fields
- `resnum::Int`: the residue number (an `AtomKey.resnum`) of the axis atoms
  `b` and `c`.
- `name::Union{Nothing,Symbol}`: `:ψ` or `:φ` for a backbone dihedral,
  `nothing` for every other rotation.
- `atoms::NTuple{4,AtomKey}`: the atoms `a`, `b`, `c`, `d` of the dihedral
  `a–b–c–d`. The rotation is about the `b–c` bond, and `d` is the atom that
  the corresponding build step places.

Side-chain rotations are left unnamed because the build tables define them
from reference atoms that need not be the IUPAC χ reference atoms: lysine's
rotation about the Cα–Cβ bond, for instance, is measured here as
C–Cα–Cβ–Cγ rather than the IUPAC N–Cα–Cβ–Cγ, so a `:χ1` label would misstate
the value. The `atoms` field defines the angle exactly.
"""
struct DihedralLabel
    resnum::Int
    name::Union{Nothing,Symbol}
    atoms::NTuple{4,AtomKey}
end

function Base.show(io::IO, label::DihedralLabel)
    print(io, "DihedralLabel(", label.resnum, ", ", repr(label.name), ", ")
    join(io, ("$(a.aname)($(a.resnum))" for a in label.atoms), "-")
    print(io, ")")
end

# `:ψ` rotates about a Cα–C bond and places the next residue's N, or the
# carboxyl terminus's OXT; `:φ` rotates about an N–Cα bond and places that
# residue's own C. Every other rotation is unnamed.
function backbonename(a::AtomKey, b::AtomKey, c::AtomKey, d::AtomKey)
    anames = (a.aname, b.aname, c.aname, d.aname)
    if anames == (:N, :CA, :C, :N) && d.resnum != c.resnum
        return :ψ
    elseif anames == (:N, :CA, :C, :OXT)
        return :ψ
    elseif anames == (:C, :N, :CA, :C) && a.resnum != b.resnum
        return :φ
    end
    return nothing
end

"""
    labels = dihedrallabels(bp::BondParametrization)

Name the rotatable dihedrals of `bp`, one `DihedralLabel` per entry of a
`dihedrals` vector and hence per column of the coordinate Jacobian. The
result has length `ndihedrals(bp)` and follows the same order as
`dihedralangles`.

Backbone dihedrals are named `:ψ` and `:φ`; side-chain dihedrals are
unnamed, and a label's `atoms` field is what defines the angle. Labels are
distinct from one another and can be used as dictionary keys.
"""
function dihedrallabels(bp::BondParametrization)
    labels = Vector{DihedralLabel}(undef, ndihedrals(bp))
    k = 0
    for step in bp.steps
        (step isa Extend && step.rotatable) || continue
        ia, ib, ic = step.predecessors
        a, b, c, d = bp.atoms[ia], bp.atoms[ib], bp.atoms[ic], bp.atoms[step.aidx]
        labels[k+=1] = DihedralLabel(c.resnum, backbonename(a, b, c, d), (a, b, c, d))
    end
    return labels
end

"""
    bp, dihedrals = bondparametrization(chain)
    bp, dihedrals = bondparametrization(T, chain)

Represent a protein `chain` by fixed bond parameters `bp` and its rotatable
`dihedrals`. These values determine its coordinates up to a rigid transform.

`dihedrals` is a vector representing the rotatable dihedral angles (in
radians) in the chain, in the order documented in [Conventions](@ref);
`ndihedrals` gives its length, `dihedrallabels` names its entries, and it is
what `dihedralangles(bp, chain)` measures. `bp` is a `BondParametrization`
object containing all other necessary information.

The first form stores `bp`'s bond lengths, bond angles, and fixed dihedrals
(and returns `dihedrals`) as `Float64`; pass a type explicitly, e.g.
`bondparametrization(Float32, chain)`, to use a different real type instead.

An `ArgumentError` reports a chain the build tables cannot describe: an
empty chain, an unrecognized residue name, a residue carrying an insertion
code, a missing backbone or side-chain atom, a chain break that leaves a
`-`/`+` build-table reference unresolvable, or a residue with atoms the
build tables never place.
"""
function bondparametrization(::Type{T}, chain::Chain) where {T<:Real}
    ress = collectresidues(chain)
    nres = length(ress)
    nres >= 1 || throw(ArgumentError("chain has no residues"))
    natoms = sum(length(res) for res in ress)
    atoms = Vector{AtomKey}(undef, natoms)
    steps = Vector{Union{Extend{T},Branch{T}}}()
    X0 = Vector{SVector{3,Float64}}(undef, natoms)   # the chain's coordinates in `atoms` order

    atomidx = Dict{String,Int}()
    resatomidxs = Vector{typeof(atomidx)}(undef, nres)
    aidx = 0
    function addatom(a::Atom)
        aidx += 1
        atoms[aidx] = AtomKey(a)
        X0[aidx] = SVector{3}(a.coords)
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
        n, cα, c = backboneatom(res, "N"), backboneatom(res, "CA"), backboneatom(res, "C")
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
            # Cα_i is fixed by the peptide bond's ω.
            push!(steps, Extend{T}((previdx[2], previdx[3], idxn), idxcα, norm(n.coords - cα.coords),
                                   bondangle(prevc, n, cα), false, ωs[i]))
            # C_i rotates about the N_i–Cα_i bond: φ_i. A ring through N–Cα
            # fixes φ (as in proline).
            rot = get(residue_phi_rotatable, resname(res), nothing)
            rot === nothing && throw(ArgumentError(unrecognized_residue_message(res)))
            push!(steps, Extend{T}((previdx[3], idxn, idxcα), idxc, norm(cα.coords - c.coords),
                                   bondangle(n, cα, c), rot, ϕs[i]))
        end
        if i == nres
            # OXT closes the carboxyl terminus; its dihedral is a final ψ.
            oxt = backboneatom(res, "OXT")
            idxoxt = addatom(oxt)
            ψ′ = dihedralangle(n, cα, c, oxt)
            push!(steps, Extend{T}((idxn, idxcα, idxc), idxoxt, norm(c.coords - oxt.coords),
                                   bondangle(cα, c, oxt), true, ψ′))
        end
        previdx = (idxn, idxcα, idxc)
        resatomidxs[i] = copy(atomidx)
    end
    # Sidechain atoms
    for (i, res) in enumerate(ress)
        rname = resname(res)
        # CYX is disulfide-bonded cysteine and therefore has no thiol hydrogen.
        rname ∈ ("CYX", "NCYX", "CCYX") && haskey(res.atoms, "HG") &&
            throw(ArgumentError("$(resdesc(res)): CYX (disulfide-bonded cysteine) must not have a thiol HG"))
        seq = get(residue_build_sequence, rname, nothing)
        seq === nothing && throw(ArgumentError(unrecognized_residue_message(res)))
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
                atomD = templateatom(res, d)
                ℓcd = norm(atomC.coords - atomD.coords)
                θbcd = bondangle(atomB, atomC, atomD)
                ϕ = dihedralangle(atomA, atomB, atomC, atomD)
                push!(steps, Extend{T}(predecessors, addatom(atomD), ℓcd, θbcd, rotatable, ϕ))
            elseif isa(step, Tuple{String,String,String,Vector{String}})
                # Branch
                a, b, c, ats = step
                i == nres && length(ats) == 1 && only(ats) == "OXT" && continue # already added
                atomA, atomB, atomC = templateatom(res, a), templateatom(res, b), templateatom(res, c)
                atomsD = [templateatom(res, at) for at in ats]
                predecessors = (atomidx[a], atomidx[b], atomidx[c])
                βs = betas(SVector{3}(atomA.coords), SVector{3}(atomB.coords), SVector{3}(atomC.coords),
                           [SVector{3}(at.coords) for at in atomsD])
                push!(steps, Branch{T}(predecessors, aidx + 1, βs))
                for at in atomsD
                    addatom(at)
                end
            else
                error("Invalid step: $step")
            end
        end
    end
    if aidx != natoms
        placed = Set(view(atoms, 1:aidx))
        detail = ""
        for res in ress
            unplaced = sort!([name for name in keys(res.atoms) if AtomKey(res[name]) ∉ placed])
            if !isempty(unplaced)
                detail = "; $(resdesc(res)) has atom(s) " * join(unplaced, ", ") *
                         " that residue_build_sequence never places"
                break
            end
        end
        throw(ArgumentError("the build sequence placed $aidx atoms but the chain has $natoms$detail"))
    end
    ndih = count(step -> step isa Extend && step.rotatable, steps)
    bp = BondParametrization{T}(atoms, nres, ℓnca, ℓcac, θncac, steps, ndih)
    return bp, convert(Vector{T}, dihedralangles(bp, X0))
end
bondparametrization(chain::Chain) = bondparametrization(Float64, chain)

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
