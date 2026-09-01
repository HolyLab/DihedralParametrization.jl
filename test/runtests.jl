using DihedralParametrization
using Aqua
using ExplicitImports
using BioStructures
using StaticArrays
using LinearAlgebra
using ForwardDiff
using Random
using Test

const testad = Base.VERSION.major == 1 && Base.VERSION.minor == 10
if testad
    using DifferentiationInterface
    using Mooncake: Mooncake
    using FiniteDifferences: central_fdm
end

# Residue `resnum`'s own φ is the dihedral of the `Extend` step that places its C.
function rotatablephi(bp, resnum)
    k = findfirst(==(DihedralParametrization.AtomData(resnum, :C)), bp.atoms)
    k === nothing && error("no C atom for residue $resnum")
    j = findfirst(s -> s isa DihedralParametrization.Extend && s.aidx == k, bp.steps)
    return bp.steps[j].rotatable
end

# Steps that place an atom other than a backbone N, Cα, C, or OXT.
issidechain(bp, step) = bp.atoms[step.aidx].aname ∉ (:N, :CA, :C, :OXT)

@testset "DihedralParametrization.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(DihedralParametrization)
    end

    @testset "ExplicitImports" begin
        # The public-ness checks fall back to `isexported` before Julia 1.11,
        # where they false-positive on `public`-but-unexported bindings.
        test_explicit_imports(DihedralParametrization;
                              all_explicit_imports_are_public   = VERSION >= v"1.11",
                              all_qualified_accesses_are_public = VERSION >= v"1.11")
    end

    @testset "Geometry" begin
        for _ = 1:5
            a, b, c, dgt = randn(SVector{3,Float64}), randn(SVector{3,Float64}), randn(SVector{3,Float64}), randn(SVector{3,Float64})
            θ = bondangle(b - c, dgt - c)
            ϕ = dihedralangle(b-a, c-b, dgt-c)
            ℓ = norm(dgt - c)
            drc = DihedralParametrization.snnerf(a, b, c, ℓ, θ, ϕ)
            @test isapprox(drc, dgt; atol=1e-8)
            βs = DihedralParametrization.betas(a, b, c, [dgt])
            drc = only(DihedralParametrization.add_to_middle(a, b, c, βs))
            @test isapprox(drc, dgt; atol=1e-8)
        end
    end

    @testset "Roundtrip" begin
        # Load a reasonably small protein with all 20 amino acids
        # (Hydrogens were added by ChimeraX)
        path = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens.cif")
        struc = read(path, MMCIFFormat)
        chain = only(only(struc))
        specializeresnames!(struc)
        bp, dihedrals = bondparametrization(chain)
        X = atomcoordinates(bp, dihedrals, chain)
        for (i, (adata, coords)) in enumerate(zip(bp.atoms, X))
            at = chain[adata.ridx][String(adata.aname)]
            @test isapprox(at.coords, coords; atol=1e-8)
        end

        # Reference coordinates as plain vectors, views, and mixed containers
        r1 = chain[1]
        n, cα, c = r1["N"].coords, r1["CA"].coords, r1["C"].coords
        @test atomcoordinates(bp, dihedrals, n, cα, c) == X
        @test atomcoordinates(bp, dihedrals, view(n, :), SVector{3}(cα), c) == X
        @test_throws DimensionMismatch atomcoordinates(bp, dihedrals, n[1:2], cα, c)
    end

    @testset "dihedralangles" begin
        path = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens.cif")
        struc = read(path, MMCIFFormat)
        specializeresnames!(struc)
        chain = only(only(struc))
        bp, dihedrals = bondparametrization(chain)

        # `bondparametrization` reports exactly what `dihedralangles` measures.
        @test dihedralangles(bp, chain) == dihedrals
        @test dihedralangles(bp, atomcoordinates(bp, dihedrals, chain)) ≈ dihedrals

        # Round trip through `atomcoordinates`, modulo 2π.
        rng = Random.Xoshiro(23)
        θ = 2π .* (rand(rng, length(dihedrals)) .- 0.5)
        θback = dihedralangles(bp, atomcoordinates(bp, θ, chain))
        @test all(x -> isapprox(x, 0; atol=1e-10), rem2pi.(θ .- θback, RoundNearest))

        # A Float32 parametrization with Float64 coordinates promotes.
        bp32, _ = bondparametrization(Float32, chain)
        X = atomcoordinates(bp, dihedrals, chain)
        @test dihedralangles(bp32, X) isa Vector{Float64}
        @test dihedralangles(bp32, chain) isa Vector{Float64}

        @test_throws DimensionMismatch dihedralangles(bp, X[1:end-1])
        @test_throws DimensionMismatch dihedralangles(bp, push!(copy(X), X[end]))

        # A chain whose atom count differs from the parametrization's.
        strucshort = read(path, MMCIFFormat)
        specializeresnames!(strucshort)
        chainshort = only(only(strucshort))
        resshort = collectresidues(chainshort)[5]
        delete!(resshort.atoms, "HB2")
        deleteat!(resshort.atom_list, findfirst(==("HB2"), resshort.atom_list))
        @test_throws "atoms but the parametrization describes" dihedralangles(bp, chainshort)

        # A chain with the right atom count but an atom the parametrization
        # cannot find: renumbering a residue moves all of its atoms.
        strucmoved = read(path, MMCIFFormat)
        specializeresnames!(strucmoved)
        chainmoved = only(only(strucmoved))
        collectresidues(chainmoved)[5].number = 10005
        @test_throws "chain has no atom N in residue 5" dihedralangles(bp, chainmoved)
    end

    @testset "dihedrallabels" begin
        path = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens.cif")
        struc = read(path, MMCIFFormat)
        specializeresnames!(struc)
        chain = only(only(struc))
        ress = collectresidues(chain)
        nres = length(ress)
        bp, dihedrals = bondparametrization(chain)
        labels = dihedrallabels(bp)

        @test length(labels) == DihedralParametrization.ndihedrals(bp)
        @test axes(labels) == axes(dihedrals)

        # The first dihedral is residue 1's psi.
        AtomData = DihedralParametrization.AtomData
        @test labels[1] == DihedralParametrization.DihedralLabel(
            1, :ψ, (AtomData(1, :N), AtomData(1, :CA), AtomData(1, :C), AtomData(2, :N)))
        @test occursin("N(1)-CA(1)-C(1)-N(2)", sprint(show, labels[1]))

        # One psi per residue: nres-1 peptide psis plus the terminal OXT.
        @test count(l -> l.name === :ψ, labels) == nres
        @test count(l -> l.name === :φ, labels) ==
              count(i -> rotatablephi(bp, resnumber(ress[i])), 2:nres)

        # Only the steps placing a backbone N, C, or OXT are named.
        @test all(l -> (l.name === nothing) == (l.atoms[4].aname ∉ (:N, :C, :OXT)), labels)

        # A ring-constrained phi contributes no label.
        prolinenums = [resnumber(r) for r in ress if resname(r) in ("PRO", "NPRO", "CPRO")]
        @test !isempty(prolinenums)
        @test !any(l -> l.name === :φ && l.resnum in prolinenums, labels)

        # Lysine's rotation about Cα–Cβ is measured from C, not from N.
        lysnum = resnumber(first(filter(r -> resname(r) in ("LYS", "NLYS", "CLYS"), ress)))
        chi1 = only(filter(l -> l.resnum == lysnum && l.atoms[4].aname === :CG, labels))
        @test chi1.name === nothing
        @test chi1.atoms == (AtomData(lysnum, :C), AtomData(lysnum, :CA),
                             AtomData(lysnum, :CB), AtomData(lysnum, :CG))

        # Labels are distinct and hashable.
        @test length(unique(labels)) == length(labels)
        bydihedral = Dict(labels .=> dihedrals)
        @test length(bydihedral) == length(labels)
        @test bydihedral[labels[7]] == dihedrals[7]

        # The labels depend only on the build sequence, not the element type.
        bp32, _ = bondparametrization(Float32, chain)
        @test dihedrallabels(bp32) == labels
    end

    @testset "Element types, inference, and corrupt parametrizations" begin
        path = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens.cif")
        struc = read(path, MMCIFFormat)
        specializeresnames!(struc)
        chain = only(only(struc))
        bp, dihedrals = bondparametrization(chain)
        X = atomcoordinates(bp, dihedrals, chain)
        r1 = first(collectresidues(chain))
        n, ca, c = SVector{3}(r1["N"].coords), SVector{3}(r1["CA"].coords), SVector{3}(r1["C"].coords)

        # A Float32 parametrization with Float64 references promotes to Float64.
        bp32, d32 = bondparametrization(Float32, chain)
        X32 = atomcoordinates(bp32, d32, chain)
        @test length(X32) == length(bp32.atoms)
        @test eltype(X32) === SVector{3,Float64}
        @test isapprox(X32, X; rtol=1e-4)

        # Dual-number references propagate into the result.
        dual(x) = SVector{3}(ForwardDiff.Dual.(x, 1.0))
        Xdual = atomcoordinates(bp, dihedrals, dual(n), dual(ca), dual(c))
        @test length(Xdual) == length(bp.atoms)
        @test eltype(Xdual) <: SVector{3,<:ForwardDiff.Dual}
        @test all(x -> isapprox(ForwardDiff.value.(x[1]), x[2]; atol=1e-8), zip(Xdual, X))

        @test @inferred(atomcoordinates(bp, dihedrals, n, ca, c)) == X
        dualdihedrals = ForwardDiff.Dual.(dihedrals, 1.0)
        @test eltype(@inferred(atomcoordinates(bp, dualdihedrals, n, ca, c))) <: SVector{3,<:ForwardDiff.Dual}

        # `dihedrals` may be any AbstractVector.
        @test atomcoordinates(bp, view(dihedrals, :), chain) == X

        # A `bp` whose steps leave an atom unplaced is rejected by both
        # traversals.
        badbp = DihedralParametrization.BondParametrization{Float64}(
            bp.atoms, bp.nres, bp.ℓnca, bp.ℓcac, bp.θncac, bp.steps[1:end-1], bp.ndihedrals)
        @test_throws "atoms but bp.atoms has" atomcoordinates(badbp, dihedrals, chain)
        @test_throws "atoms but bp.atoms has" jacobianplan(badbp)
    end

    @testset "buildchain" begin
        path = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens.cif")
        struc = read(path, MMCIFFormat)
        specializeresnames!(struc)
        chain = only(only(struc))
        bp, dihedrals = bondparametrization(chain)

        X = atomcoordinates(bp, dihedrals, chain)
        rebuilt = buildchain(chain, bp, X)
        for adata in bp.atoms
            a0 = chain[adata.ridx][String(adata.aname)]
            a1 = rebuilt[adata.ridx][String(adata.aname)]
            @test isapprox(a0.coords, a1.coords; atol=1e-8)
        end
        _, dihedralsround = bondparametrization(rebuilt)
        @test isapprox(dihedralsround, dihedrals; atol=1e-8)

        # Perturbed dihedrals: rebuilding and re-encoding must recover the
        # same rotatable dihedrals modulo 2π.
        rng = Random.Xoshiro(7)
        θpert = dihedrals .+ 0.3 .* randn(rng, length(dihedrals))
        Xpert = atomcoordinates(bp, θpert, chain)
        pertchain = buildchain(chain, bp, Xpert)
        _, dihedralspert = bondparametrization(pertchain)
        wrap(x) = mod(x + π, 2π) - π
        @test all(k -> isapprox(wrap(dihedralspert[k] - θpert[k]), 0; atol=1e-8), eachindex(dihedrals))
    end

    @testset "2QMT (experimental structure)" begin
        # Provenance: test/data/README.md.
        path = joinpath(@__DIR__, "data", "2qmt_H.cif")
        struc = read(path, MMCIFFormat)
        specializeresnames!(struc)
        chain = only(only(struc))
        bp, dihedrals = bondparametrization(chain)

        X = atomcoordinates(bp, dihedrals, chain)
        for (adata, coords) in zip(bp.atoms, X)
            at = chain[adata.ridx][String(adata.aname)]
            @test isapprox(at.coords, coords; atol=1e-8)
        end

        rebuilt = buildchain(chain, bp, X)
        for adata in bp.atoms
            a0 = chain[adata.ridx][String(adata.aname)]
            a1 = rebuilt[adata.ridx][String(adata.aname)]
            @test isapprox(a0.coords, a1.coords; atol=1e-8)
        end
    end

    @testset "Disulfide (CYX)" begin
        # Residues 37 and 69 are renamed CYX and have no HG. The fixture does
        # not model the S-S bond because this test only exercises templates.
        path = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens_CYX.cif")
        struc = read(path, MMCIFFormat)
        chain = only(only(struc))
        specializeresnames!(struc)
        ress = collectresidues(chain)
        @test count(r -> resname(r) == "CYX", ress) == 2

        bp, dihedrals = bondparametrization(chain)
        X = atomcoordinates(bp, dihedrals, chain)
        for (adata, coords) in zip(bp.atoms, X)
            at = chain[adata.ridx][String(adata.aname)]
            @test isapprox(at.coords, coords; atol=1e-8)
        end
        cyxridxs = Set(resnumber(r) for r in ress if resname(r) == "CYX")
        @test !any(a -> a.ridx in cyxridxs && a.aname == :HG, bp.atoms)

        # CYX must not retain its thiol hydrogen.
        pathfree = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens.cif")
        strucfree = read(pathfree, MMCIFFormat)
        chainfree = only(only(strucfree))
        specializeresnames!(strucfree)
        ressfree = collectresidues(chainfree)
        cysfree = only(filter(r -> resnumber(r) == 37, ressfree))
        @test resname(cysfree) == "CYS"
        cysfree.name = "CYX"   # CYX with HG still attached
        @test_throws "CYX (disulfide-bonded cysteine) must not have a thiol HG" bondparametrization(chainfree)

        # Missing template atoms report residue context.
        cysfree.name = "CYS"
        delete!(cysfree.atoms, "HG")
        @test_throws "cannot find atom \"HG\"" bondparametrization(chainfree)
    end

    @testset "Bond lengths and angles are invariant under every dihedral" begin
        # Changing dihedrals must preserve all bond lengths and angles,
        # including those spanning residue boundaries.
        path = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens.cif")
        struc = read(path, MMCIFFormat)
        specializeresnames!(struc)
        chain = only(only(struc))
        ress = collectresidues(chain)
        nres = length(ress)
        bp, dihedrals = bondparametrization(chain)

        X0 = atomcoordinates(bp, dihedrals, chain)
        for (adata, coords) in zip(bp.atoms, X0)
            at = chain[adata.ridx][String(adata.aname)]
            @test isapprox(at.coords, coords; atol=1e-8)
        end

        natoms = length(X0)
        anames = [String(a.aname) for a in bp.atoms]
        ridx = [a.ridx for a in bp.atoms]

        # Derive connectivity from residue templates, not atom distances.
        atomrow = Dict((a.ridx, String(a.aname)) => k for (k, a) in pairs(bp.atoms))
        bonds = Tuple{Int,Int}[]
        for res in ress
            rd = BioStructures.residuedata[resname(res)]
            for (n1, n2) in rd.bonds
                (haskey(res.atoms, n1) && haskey(res.atoms, n2)) || continue
                push!(bonds, (atomrow[(resnumber(res), n1)], atomrow[(resnumber(res), n2)]))
            end
        end
        for i = 1:nres-1     # peptide bonds
            push!(bonds, (atomrow[(resnumber(ress[i]), "C")], atomrow[(resnumber(ress[i+1]), "N")]))
        end
        nbrs = [Int[] for _ = 1:natoms]
        for (i, j) in bonds
            push!(nbrs[i], j)
            push!(nbrs[j], i)
        end
        angles = Tuple{Int,Int,Int}[]
        for b = 1:natoms, ii in eachindex(nbrs[b]), jj = ii+1:length(nbrs[b])
            push!(angles, (nbrs[b][ii], b, nbrs[b][jj]))
        end

        # Confirm that cross-residue bonds and angles are covered.
        crossbond(i, j) = ridx[i] != ridx[j]
        crosstriple(a, b, c) = ridx[a] != ridx[b] || ridx[b] != ridx[c]
        namematch(i, j, n1, n2) = (anames[i] == n1 && anames[j] == n2) || (anames[i] == n2 && anames[j] == n1)
        nCN = count(((i, j),) -> crossbond(i, j) && namematch(i, j, "C", "N"), bonds)
        nONC = count(((a, b, c),) -> crosstriple(a, b, c) && anames[b] == "C" && namematch(a, c, "O", "N"), angles)
        nCNH = count(((a, b, c),) -> crosstriple(a, b, c) && anames[b] == "N" && namematch(a, c, "C", "H"), angles)
        nCNCA = count(((a, b, c),) -> crosstriple(a, b, c) && anames[b] == "N" && namematch(a, c, "C", "CA"), angles)
        @test nCN == nres - 1
        @test nONC == nres - 1
        @test nCNCA == nres - 1
        # Proline's ring nitrogen has no amide H.
        @test nCNH == count(i -> haskey(ress[i].atoms, "H"), 2:nres)

        # Bond lengths and angles must have zero derivatives with respect to
        # every dihedral.
        function internal_coords(dih)
            X = atomcoordinates(bp, dih, chain)
            ls = [norm(X[i] - X[j]) for (i, j) in bonds]
            as = [bondangle(X[a] - X[b], X[c] - X[b]) for (a, b, c) in angles]
            return vcat(ls, as)
        end
        J = ForwardDiff.jacobian(internal_coords, dihedrals)
        nbonds = length(bonds)
        tol = 1e-8
        sensitive(k) = maximum(abs, @view J[k, :]) > tol

        badlengths = count(sensitive, 1:nbonds)
        badangles = count(sensitive, nbonds+1:size(J, 1))
        @test badlengths == 0
        @test badangles == 0
    end

    @testset "Degrees of freedom" begin
        # Proline's ring fixes chi1, chi2, and phi.
        path = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens.cif")
        struc = read(path, MMCIFFormat)
        specializeresnames!(struc)
        chain = only(only(struc))
        ress = collectresidues(chain)
        nres = length(ress)
        nproline = count(r -> resname(r) in ("PRO", "NPRO", "CPRO"), ress)
        @test nproline == 8

        bp, dihedrals = bondparametrization(chain)

        # Backbone: psi for every residue but the last, phi for every residue
        # but the first and any proline, plus one final dihedral to place OXT.
        nphi_free = count(i -> rotatablephi(bp, resnumber(ress[i])), 2:nres)
        @test nphi_free == nres - 1 - nproline
        nsidechain_free = count(bp.steps) do step
            step isa DihedralParametrization.Extend && step.rotatable && issidechain(bp, step)
        end
        @test length(dihedrals) == (nres - 1) + nphi_free + 1 + nsidechain_free

        # No proline contributes a phi or a chi to `dihedrals`.
        prolineresidx = [i for i = 1:nres if resname(ress[i]) in ("PRO", "NPRO", "CPRO")]
        for i in prolineresidx
            i > 1 && @test !rotatablephi(bp, resnumber(ress[i]))
            @test !any(bp.steps) do step
                step isa DihedralParametrization.Extend && step.rotatable && issidechain(bp, step) &&
                    bp.atoms[step.aidx].ridx == resnumber(ress[i])
            end
        end
    end

    @testset "Analytic Jacobian" begin
        flatten(X) = reduce(vcat, Vector.(X))

        path = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens.cif")
        struc = read(path, MMCIFFormat)
        specializeresnames!(struc)
        chain = only(only(struc))
        ress = collectresidues(chain)
        nres = length(ress)
        bp, dihedrals = bondparametrization(chain)
        n, cα, c = SVector{3}(ress[1]["N"].coords), SVector{3}(ress[1]["CA"].coords), SVector{3}(ress[1]["C"].coords)

        plan = jacobianplan(bp)
        @test plan.ndih == length(dihedrals)

        # Compare with ForwardDiff at two configurations.
        X = atomcoordinates(bp, dihedrals, n, cα, c)
        J = coordinatejacobian(plan, X)
        Jad = ForwardDiff.jacobian(θ -> flatten(atomcoordinates(bp, collect(θ), n, cα, c)), dihedrals)
        @test isapprox(J, Jad; rtol=1e-10)

        rng = Random.Xoshiro(42)
        θpert = dihedrals .+ 0.3 .* randn(rng, length(dihedrals))
        Xpert = atomcoordinates(bp, θpert, n, cα, c)
        Jpert = coordinatejacobian(plan, Xpert)
        Jadpert = ForwardDiff.jacobian(θ -> flatten(atomcoordinates(bp, collect(θ), n, cα, c)), θpert)
        @test isapprox(Jpert, Jadpert; rtol=1e-10)

        Jsmall = zeros(size(Jpert) .- 1)
        @test_throws DimensionMismatch coordinatejacobian!(Jsmall, plan, Xpert)

        # Compare products with the explicit Jacobian.
        for _ = 1:3
            w = [randn(rng, SVector{3,Float64}) for _ = 1:plan.natoms]
            g = vjp!(zeros(plan.ndih), plan, Xpert, w)
            @test isapprox(g, Jpert' * flatten(w); rtol=1e-12)

            v = randn(rng, plan.ndih)
            δx = jvp!([zero(SVector{3,Float64}) for _ = 1:plan.natoms], plan, Xpert, v)
            @test isapprox(flatten(δx), Jpert * v; rtol=1e-12)

            # The allocating forms agree with the in-place ones.
            @test vjp(plan, Xpert, w) == g
            @test jvp(plan, Xpert, v) == δx

            # `reinterpret` relates the coordinate list to the rows of `J`.
            @test isapprox(Jpert' * reinterpret(Float64, w), vjp(plan, Xpert, w); rtol=1e-12)
            @test isapprox(reinterpret(Float64, jvp(plan, Xpert, v)), Jpert * v; rtol=1e-12)
        end

        # Element types of the allocating forms follow promotion.
        let X32 = [SVector{3,Float32}(x) for x in Xpert],
            w32 = [randn(rng, SVector{3,Float32}) for _ = 1:plan.natoms],
            v32 = randn(rng, Float32, plan.ndih),
            w64 = [randn(rng, SVector{3,Float64}) for _ = 1:plan.natoms],
            v64 = randn(rng, plan.ndih)

            @test eltype(vjp(plan, X32, w32)) === Float32
            @test eltype(vjp(plan, X32, w64)) === Float64
            @test eltype(eltype(jvp(plan, X32, v32))) === Float32
            @test eltype(eltype(jvp(plan, X32, v64))) === Float64
        end

        # The first three atoms (N, Cα, C of residue 1) never move.
        @test all(iszero, @view J[1:9, :])

        # Exercise fixed proline phi.
        prolineidx = findfirst(i -> resname(ress[i]) in ("PRO", "NPRO", "CPRO") && i > 1, 1:nres)
        @test prolineidx !== nothing
        @test !rotatablephi(bp, resnumber(ress[prolineidx]))

        # Map backbone dihedrals to columns.
        idx = 0
        psicol = Dict{Int,Int}()
        phicol = Dict{Int,Int}()
        for i = 1:nres-1
            psicol[i] = (idx += 1)
            rotatablephi(bp, resnumber(ress[i+1])) && (phicol[i+1] = (idx += 1))
        end

        # A carbonyl O with a +N reference moves under psi_i.
        atomrow = Dict((a.ridx, String(a.aname)) => k for (k, a) in pairs(bp.atoms))
        i = findfirst(i -> haskey(atomrow, (i, "O")) && haskey(psicol, i), 1:nres-1)
        Orows = 3 * (atomrow[(i, "O")] - 1) + 1 : 3 * atomrow[(i, "O")]
        @test any(!iszero, J[Orows, psicol[i]])

        # An amide H with a -C reference lies on the phi axis.
        j = findfirst(j -> haskey(atomrow, (j, "H")) && haskey(phicol, j), 2:nres)
        j = j === nothing ? nothing : j + 1
        @test j !== nothing
        Hrows = 3 * (atomrow[(j, "H")] - 1) + 1 : 3 * atomrow[(j, "H")]
        @test all(iszero, J[Hrows, phicol[j]])

        # Independent finite-difference check.
        function fdcolumn(k; h=1e-6)
            θp = copy(dihedrals); θp[k] += h
            θm = copy(dihedrals); θm[k] -= h
            return (flatten(atomcoordinates(bp, θp, n, cα, c)) .- flatten(atomcoordinates(bp, θm, n, cα, c))) ./ 2h
        end
        for k in (1, cld(plan.ndih, 2), plan.ndih)
            @test isapprox(fdcolumn(k), J[:, k]; rtol=1e-6, atol=1e-6)
        end

        # Exercise a disulfide.
        pathcyx = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens_CYX.cif")
        struccyx = read(pathcyx, MMCIFFormat)
        specializeresnames!(struccyx)
        chaincyx = only(only(struccyx))
        bpcyx, dihedralscyx = bondparametrization(chaincyx)
        plancyx = jacobianplan(bpcyx)
        @test plancyx.ndih == length(dihedralscyx)
    end

    @testset "Analytic weighted Hessian" begin
        flatten(X) = reduce(vcat, Vector.(X))

        path = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens.cif")
        struc = read(path, MMCIFFormat)
        specializeresnames!(struc)
        chain = only(only(struc))
        ress = collectresidues(chain)
        n, cα, c = SVector{3}(ress[1]["N"].coords), SVector{3}(ress[1]["CA"].coords), SVector{3}(ress[1]["C"].coords)
        bp, dihedrals = bondparametrization(chain)
        plan = jacobianplan(bp)

        rng = Random.Xoshiro(11)
        w = [randn(rng, SVector{3,Float64}) for _ = 1:plan.natoms]
        wf = flatten(w)
        objective(bp, n, cα, c, w, θ) = dot(wf, flatten(atomcoordinates(bp, collect(θ), n, cα, c)))

        # Native configuration.
        X = atomcoordinates(bp, dihedrals, n, cα, c)
        S = weightedhessian(plan, X, w)
        @test issymmetric(S)
        Had = ForwardDiff.hessian(θ -> objective(bp, n, cα, c, w, θ), dihedrals)
        devnative = maximum(abs, S .- Had)
        @test isapprox(S, Had; rtol=1e-10)

        # Perturbed configuration.
        θpert = dihedrals .+ 0.3 .* randn(rng, length(dihedrals))
        Xpert = atomcoordinates(bp, θpert, n, cα, c)
        Spert = weightedhessian(plan, Xpert, w)
        @test issymmetric(Spert)
        Hadpert = ForwardDiff.hessian(θ -> objective(bp, n, cα, c, w, θ), θpert)
        devpert = maximum(abs, Spert .- Hadpert)
        @test isapprox(Spert, Hadpert; rtol=1e-10)

        # weightedhessian! in-place matches weightedhessian.
        Sfilled = zeros(plan.ndih, plan.ndih)
        weightedhessian!(Sfilled, plan, Xpert, w)
        @test Sfilled == Spert

        # A second, independently seeded random w.
        rng2 = Random.Xoshiro(97)
        w2 = [randn(rng2, SVector{3,Float64}) for _ = 1:plan.natoms]
        w2f = flatten(w2)
        S2 = weightedhessian(plan, Xpert, w2)
        @test issymmetric(S2)
        Had2 = ForwardDiff.hessian(θ -> dot(w2f, flatten(atomcoordinates(bp, collect(θ), n, cα, c))), θpert)
        @test isapprox(S2, Had2; rtol=1e-10)

        # Dimension mismatches.
        Ssmall = zeros(plan.ndih - 1, plan.ndih - 1)
        @test_throws DimensionMismatch weightedhessian!(Ssmall, plan, Xpert, w)
        Xshort = Xpert[1:end-1]
        @test_throws DimensionMismatch weightedhessian!(Sfilled, plan, Xshort, w)
        wshort = w[1:end-1]
        @test_throws DimensionMismatch weightedhessian!(Sfilled, plan, Xpert, wshort)

        # Exercise a disulfide.
        pathcyx = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens_CYX.cif")
        struccyx = read(pathcyx, MMCIFFormat)
        specializeresnames!(struccyx)
        chaincyx = only(only(struccyx))
        bpcyx, dihedralscyx = bondparametrization(chaincyx)
        plancyx = jacobianplan(bpcyx)
        rescyx = collectresidues(chaincyx)
        ncyx = SVector{3}(rescyx[1]["N"].coords)
        cαcyx = SVector{3}(rescyx[1]["CA"].coords)
        ccyx = SVector{3}(rescyx[1]["C"].coords)
        Xcyx = atomcoordinates(bpcyx, dihedralscyx, ncyx, cαcyx, ccyx)
        wcyx = [randn(rng, SVector{3,Float64}) for _ = 1:plancyx.natoms]
        wcyxf = flatten(wcyx)
        Scyx = weightedhessian(plancyx, Xcyx, wcyx)
        @test issymmetric(Scyx)
        Hadcyx = ForwardDiff.hessian(θ -> dot(wcyxf, flatten(atomcoordinates(bpcyx, collect(θ), ncyx, cαcyx, ccyx))), dihedralscyx)
        @test isapprox(Scyx, Hadcyx; rtol=1e-10)
    end

    @testset "Error handling" begin
        path = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens.cif")
        struc = read(path, MMCIFFormat)
        specializeresnames!(struc)
        chain = only(only(struc))
        bp, dihedrals = bondparametrization(chain)

        @test_throws "rotatable dihedrals" atomcoordinates(bp, dihedrals[1:end-1], chain)
        @test_throws "rotatable dihedrals" atomcoordinates(bp, vcat(dihedrals, 0.0), chain)

        # An un-renamed HIS trips the phi-rotatability lookup, which runs
        # while the residue's own backbone is encoded.
        strucraw = read(path, MMCIFFormat)
        chainraw = only(only(strucraw))
        @test_throws "histidine must be disambiguated as HID/HIE/HIP" bondparametrization(chainraw)

        # Residue 1's name is first checked at its build-sequence lookup.
        struc1 = read(path, MMCIFFormat)
        specializeresnames!(struc1)
        chain1 = only(only(struc1))
        ress1 = collectresidues(chain1)
        ress1[1].name = "HIS"
        @test_throws "residue 1 (HIS): unrecognized residue name" bondparametrization(chain1)

        plan = jacobianplan(bp)
        X = atomcoordinates(bp, dihedrals, chain)
        w = [zero(SVector{3,Float64}) for _ = 1:plan.natoms]
        v = zeros(plan.ndih)

        @test_throws DimensionMismatch vjp!(zeros(plan.ndih - 1), plan, X, w)
        @test_throws DimensionMismatch vjp!(zeros(plan.ndih), plan, X, w[1:end-1])
        @test_throws DimensionMismatch vjp!(zeros(plan.ndih), plan, X[1:end-1], w)

        @test_throws DimensionMismatch jvp!([zero(SVector{3,Float64}) for _ = 1:plan.natoms-1], plan, X, v)
        @test_throws DimensionMismatch jvp!([zero(SVector{3,Float64}) for _ = 1:plan.natoms], plan, X, v[1:end-1])
        @test_throws DimensionMismatch jvp!([zero(SVector{3,Float64}) for _ = 1:plan.natoms], plan, X[1:end-1], v)
    end

    if testad
        @testset "Differentiability" begin
            function makef(chain)  # Keep inference check self-contained.
                bp, dihedrals = bondparametrization(chain)
                return function(dih)
                    return atomcoordinates(bp, dih, chain)
                end, length(dihedrals)
            end

            path = joinpath(@__DIR__, "data", "AF-M3YHX5-F1-model_v4_hydrogens.cif")
            struc = read(path, MMCIFFormat)
            specializeresnames!(struc)
            A = only(only(struc))

            makeX, n = makef(A)
            dih = 2π * (rand(n) .- 0.5)
            X = makeX(dih)
            backend = DifferentiationInterface.AutoMooncakeForward(; config=nothing)
            J = jacobian(makeX, backend, dih)
            Jmc = [SVector(j.fields.data) for j in J]
            backend = AutoFiniteDifferences(; fdm=central_fdm(5, 1))
            J = jacobian(makeX, backend, dih)
            Jfd = reinterpret(SVector{3,Float64}, J)
            @test Jmc ≈ Jfd atol=1e-6
        end
    end
end
