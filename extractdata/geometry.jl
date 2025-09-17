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

# The remaining cases do not have any flexibility.
# We still need to know how to add them to the growing chain
# Parse the residue graphs to determine the bonding patterns
function parse_residue_graph(resname::AbstractString)
    rd = BioStructures.residuedata[resname]
    g = SimpleGraph()
    mg = MetaGraph(g, Symbol, Int, Nothing)
    idx = 0
    # Add vertices, sorting atoms into rotatable vs. inflexible sets (the rotatable will be paired with an entry in `dihedrals`)
    # The backbone atoms are always placed first
    placed, inflexible = Set{Symbol}([:N, :CA, :C]), Set{Symbol}()
    for (atomname, _) in rd.atoms
        sym = Symbol(atomname)
        mg[sym] = idx += 1
        sym ∉ placed && push!(inflexible, sym)
    end
    lookup = resname
    if length(lookup) == 4 && lookup[1] ∈ ('N', 'C')
        lookup = lookup[2:end]
    end
    lookup = lookup ∈ ("HID", "HIE", "HIP") ? "HIS" : lookup   # in chitables, histidine is only cataloged once, as "HIS"
    for ct in chitables
        quad = get(ct, lookup, nothing)
        quad === nothing && break
        sym = Symbol(quad[4])
        delete!(inflexible, sym)
    end
    if length(resname) == 4
        if resname[1] == 'N'
            delete!(inflexible, :H2)
            delete!(inflexible, :H3)
        elseif resname[1] == 'C'
            delete!(inflexible, :OXT)
        end
    end
    # Add edges
    for bond in rd.bonds
        add_edge!(g, mg[Symbol(bond[1])], mg[Symbol(bond[2])])
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
            if a === :HA
                if c === :N && d === :H
                    additions[i] = (:C, b, c, d, tf)
                elseif c === :C && d === :O
                    additions[i] = (:N, b, c, d, tf)
                end
            end
        end
    end
    return additions
end

const generate_for_aas = Set(["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"])
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
    println(io, "const residue_build_sequence = Dict(")
    for resname in sort(collect(keys(BioStructures.residuedata)))
        resname ∉ generate_for_aas && continue
        additions = parse_residue_graph(resname)
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
end

# Sanity checks
include(joinpath(dirname(@__DIR__), "src", "tables.jl"))
for (key, list) in residue_build_sequence
    nflex = sum(item -> last(item) === true, list)
    key0 = key
    offset = 0
    if length(key) == 4 && key[1] ∈ ('N', 'C')
        key = key[2:end]
        offset = key0 == "NPRO" ? 0 : 1   # PRO has no rotatable dihedrals at N because of the ring
    end
    if key ∈ ("HID", "HIE", "HIP")
        key = "HIS"
    end
    @assert nflex == something(findlast(ct -> haskey(ct, key), chitables), 0) + offset "Mismatch in number of flexible dihedrals for residue $key0"
end
