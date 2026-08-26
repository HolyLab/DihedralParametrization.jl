using BioStructures
using Graphs
using MetaGraphsNext
using LinearAlgebra

const chitables = [copy(d) for d in BioStructures.chitables]

function add_to_tables!(ct, (key, val))
    i = findfirst(d -> !haskey(d, key), ct)
    push!(ct[i], key => val)
end

# Add χ angles for added hydrogens
add_to_tables!(chitables, ("LYS", ("CD", "CE", "NZ", "HZ1")))
add_to_tables!(chitables, ("SER", ("CA", "CB", "OG", "HG")))
add_to_tables!(chitables, ("THR", ("CA", "CB", "OG1", "HG1")))
add_to_tables!(chitables, ("THR", ("CA", "CB", "CG2", "HG21")))
add_to_tables!(chitables, ("ASN", ("CB", "CG", "ND2", "HD21")))
add_to_tables!(chitables, ("GLN", ("CG", "CD", "NE2", "HE21")))
add_to_tables!(chitables, ("CYS", ("CA", "CB", "SG", "HG")))
add_to_tables!(chitables, ("ALA", ("N", "CA", "CB", "HB1")))
add_to_tables!(chitables, ("VAL", ("CA", "CB", "CG1", "HG11")))
add_to_tables!(chitables, ("VAL", ("CA", "CB", "CG2", "HG21")))
add_to_tables!(chitables, ("ILE", ("CA", "CB", "CG2", "HG21")))
add_to_tables!(chitables, ("ILE", ("CB", "CG1", "CD1", "HD11")))
add_to_tables!(chitables, ("LEU", ("CB", "CG", "CD1", "HD11")))
add_to_tables!(chitables, ("LEU", ("CB", "CG", "CD2", "HD21")))
add_to_tables!(chitables, ("MET", ("CG", "SD", "CE", "HE1")))
add_to_tables!(chitables, ("TYR", ("CE1", "CZ", "OH", "HH")))

const Branched = Tuple{Symbol,Symbol,Symbol,Vector{Symbol}}
const Extended = Tuple{Symbol,Symbol,Symbol,Symbol,Bool}

function get3(list, k::Tuple{Symbol,Symbol,Symbol})
    for s in list
        s[1] == k[1] && s[2] == k[2] && s[3] == k[3] && return s[4]
    end
    return nothing
end

"""
    g, mg, isbridge = residuegraph(resname::AbstractString)

Build the bonded-atom graph for `resname` from `BioStructures.residuedata`.
`isbridge(u::Symbol, v::Symbol)` reports whether the bond between atoms `u`
and `v` is a bridge, i.e. not part of any ring.
"""
function residuegraph(resname::AbstractString)
    rd = BioStructures.residuedata[resname]
    g = SimpleGraph()
    mg = MetaGraph(g, Symbol, Int, Nothing)
    idx = 0
    for (atomname, _) in rd.atoms
        mg[Symbol(atomname)] = idx += 1
    end
    for bond in rd.bonds
        add_edge!(g, mg[Symbol(bond[1])], mg[Symbol(bond[2])])
    end
    bridgeset = Set((min(src(e), dst(e)), max(src(e), dst(e))) for e in bridges(g))
    isbridge(u::Symbol, v::Symbol) = (min(mg[u], mg[v]), max(mg[u], mg[v])) in bridgeset
    return g, mg, isbridge
end

# Returns `(additions, phirotatable)`: `additions` is the sidechain build
# sequence, and `phirotatable` reports whether this residue's own N-CA bond is
# free to rotate (a bridge in the residue graph) rather than fixed by a ring.
function parse_residue_graph(resname::AbstractString)
    rd = BioStructures.residuedata[resname]
    g, mg, isbridge = residuegraph(resname)
    # Rotatable atoms receive entries in `dihedrals`; backbone atoms come first.
    placed, inflexible = Set{Symbol}([:N, :CA, :C]), Set{Symbol}()
    for (atomname, _) in rd.atoms
        sym = Symbol(atomname)
        sym ∉ placed && push!(inflexible, sym)
    end
    # A chi is free when its central bond is a bridge in the residue graph.
    lookup = resname
    if length(lookup) == 4 && lookup[1] ∈ ('N', 'C')
        lookup = lookup[2:end]
    end
    lookup = lookup ∈ ("HID", "HIE", "HIP") ? "HIS" : lookup
    lookup = lookup == "CYX" ? "CYS" : lookup
    for ct in chitables
        quad = get(ct, lookup, nothing)
        quad === nothing && break
        b, c, sym = Symbol(quad[2]), Symbol(quad[3]), Symbol(quad[4])
        isbridge(b, c) && delete!(inflexible, sym)
    end
    if length(resname) == 4
        if resname[1] == 'N'
            delete!(inflexible, :H2)
            delete!(inflexible, :H3)
        elseif resname[1] == 'C'
            delete!(inflexible, :OXT)
        end
    end

    # When ordering the neighbors, first deal with the rotatable ones, because the remainder may be added via a branch
    function byrot(a::Symbol)
        s = string(a)
        c = first(s)
        return (a ∈ inflexible, findfirst(==(c), "CNOSH"), length(s) >= 2 ? s[2] : '0')
    end
    function byplaced(a::Symbol)
        s = string(a)
        c = first(s)
        return (a ∈ placed, findfirst(==(c), "CNOSH"), length(s) >= 2 ? s[2] : '0')
    end

    additions = Union{Branched, Extended}[]
    lbls = collect(labels(mg))
    while length(placed) < nv(g)
        changed = false
        # First try to add branches (need two placed neighbors)
        for mid in placed
            v = mg[mid]
            ns = lbls[neighbors(g, v)]
            nsp = sort!(ns ∩ placed; by=byplaced)
            if length(nsp) >= 2
                for nbr in setdiff(ns, placed)
                    k = (nsp[1], mid, nsp[2])
                    list = get3(additions, k)
                    if list !== nothing
                        push!(list, nbr)
                    else
                        push!(additions, (k..., [nbr]))
                    end
                    push!(placed, nbr)
                    changed = true
                end
            end
        end
        # If we failed to add any branches, try to add extends (need a -> b -> c all placed)
        if length(placed) < nv(g) && !changed
            for d in sort!(setdiff(lbls, placed); by=byrot)
                nsd = lbls[neighbors(g, mg[d])] ∩ placed
                isempty(nsd) && continue
                c = first(nsd)
                nsc = lbls[neighbors(g, mg[c])] ∩ placed
                isempty(nsc) && continue
                b = only(nsc)   # `only` because otherwise we would have added a branch above
                nsb = sort!(lbls[neighbors(g, mg[b])]; by=byplaced)
                a = nothing
                for aa in nsb
                    aa == c && continue
                    a = aa
                    break
                end
                a === nothing && continue
                push!(additions, (a, b, c, d, d ∉ inflexible))
                push!(placed, d)
                changed = true
                break
            end
        end
    end
    @assert length(placed) == nv(g) || changed "Could not fully parse residue $resname; stuck with placed=$(collect(placed)), inflexible=$(collect(inflexible))"
    for (i, quad) in enumerate(additions)
        if isa(quad, Branched)
            sort!(quad[4])
        else
            a, b, c, d, tf = quad
            # Fixed dihedrals about N-CA or CA-C need a reference in the
            # adjacent residue; a same-residue reference would co-rotate.
            # Terminal variants retain their same-residue reference.
            if !tf && b === :CA && c === :N && !(length(resname) == 4 && resname[1] == 'N')
                additions[i] = (Symbol("-C"), b, c, d, tf)
            elseif !tf && b === :CA && c === :C && !(length(resname) == 4 && resname[1] == 'C')
                additions[i] = (Symbol("+N"), b, c, d, tf)
            end
        end
    end
    return additions, isbridge(:CA, :N)
end

"""
    nfreechis(resname::AbstractString, lookup::AbstractString=resname)

Count chi entries whose central bond is a bridge in `resname`'s residue graph.
Ignore entries whose target atom is absent.
"""
function nfreechis(resname::AbstractString, lookup::AbstractString=resname)
    _, mg, isbridge = residuegraph(resname)
    n = 0
    for ct in chitables
        quad = get(ct, lookup, nothing)
        quad === nothing && break
        Symbol(quad[4]) ∈ labels(mg) && isbridge(Symbol(quad[2]), Symbol(quad[3])) && (n += 1)
    end
    return n
end

const generate_for_aas = Set(["ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLU", "GLN", "GLY", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"])
for aa in collect(generate_for_aas)
    push!(generate_for_aas, "N"*aa)
    push!(generate_for_aas, "C"*aa)
end

qt(s::Symbol) = "\"$(s)\""

open(joinpath(dirname(@__DIR__), "src", "tables.jl"), "w") do io
    println(io, "# This file is auto-generated by extractdata/geometry.jl; do not edit directly.")
    println(io, "\n# Sequence of construction of atoms in residues")
    println(io, "# The pattern (a, b, c, [d1, d2, ...]) indicates that atoms d1, d2, ... are attached to b")
    println(io, "# The pattern (a, b, c, d, tf) indicates that atom d is attached to c; tf indicates whether the dihedral")
    println(io, "#   formed by a-b-c-d is rotatable (true) or fixed (false)")
    println(io, "# A reference atom prefixed with \"-\" or \"+\" (e.g. \"-C\", \"+N\") names that atom in the previous or")
    println(io, "#   next residue rather than the current one, following the CHARMM convention.")
    println(io, "const residue_build_sequence = Dict(")
    phirotatable = Dict{String,Bool}()
    for resname in sort(collect(keys(BioStructures.residuedata)))
        resname ∉ generate_for_aas && continue
        additions, phirotatable[resname] = parse_residue_graph(resname)
        println(io, "    \"$resname\" => [")
        first = true
        for quad in additions
            if isa(quad, Branched)
                a, b, c, ds = quad
                println(io, "            (\"$a\", \"$b\", \"$c\", [", join(qt.(ds), ", "), "]),")
            else
                a, b, c, d, tf = quad
                println(io, "            (\"$a\", \"$b\", \"$c\", \"$d\", $(tf)),")
            end
            first = false
        end
        println(io, "        ],")
    end
    println(io, ")")
    println(io, "\n# Whether a residue's own phi dihedral (about its N-CA bond) is a free")
    println(io, "# parameter. It is fixed when N-CA lies on a ring, e.g. proline's; the")
    println(io, "# ring's fixed bond lengths and angles then leave phi no freedom to vary.")
    println(io, "const residue_phi_rotatable = Dict(")
    for resname in sort(collect(keys(phirotatable)))
        println(io, "    \"$resname\" => $(phirotatable[resname]),")
    end
    println(io, ")")
end

# Sanity checks
include(joinpath(dirname(@__DIR__), "src", "tables.jl"))
for (key, list) in residue_build_sequence
    nflex = sum(item -> last(item) === true, list)
    key0 = key
    offset = 0
    if length(key) == 4 && key[1] ∈ ('N', 'C')
        key = key[2:end]
        offset = 1   # the added N-terminal H2 spin, or C-terminal OXT spin, is rotatable
    end
    lookup = key ∈ ("HID", "HIE", "HIP") ? "HIS" : key
    lookup = lookup == "CYX" ? "CYS" : lookup
    @assert nflex == nfreechis(key, lookup) + offset "Mismatch in number of flexible dihedrals for residue $key0"
end
