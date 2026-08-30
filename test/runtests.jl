using DihedralParametrization
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

@testset "DihedralParametrization.jl" begin
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
        nphi_free = count(i -> bp.phirotatable[i], 2:nres)
        @test nphi_free == nres - 1 - nproline
        nsidechain_free = sum(r -> count(step -> step isa DihedralParametrization.Extend && step.rotatable, r.steps), bp.residues)
        @test length(dihedrals) == (nres - 1) + nphi_free + 1 + nsidechain_free

        # No proline contributes a phi or a chi to `dihedrals`.
        prolineresidx = [i for i = 1:nres if resname(ress[i]) in ("PRO", "NPRO", "CPRO")]
        for i in prolineresidx
            i > 1 && @test !bp.phirotatable[i]
            @test !any(step -> step isa DihedralParametrization.Extend && step.rotatable, bp.residues[i].steps)
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
            g = jtv!(zeros(plan.ndih), plan, Xpert, w)
            @test isapprox(g, Jpert' * flatten(w); rtol=1e-12)

            v = randn(rng, plan.ndih)
            δx = jvp!([zero(SVector{3,Float64}) for _ = 1:plan.natoms], plan, Xpert, v)
            @test isapprox(flatten(δx), Jpert * v; rtol=1e-12)
        end

        # The first three atoms (N, Cα, C of residue 1) never move.
        @test all(iszero, @view J[1:9, :])

        # Exercise fixed proline phi.
        prolineidx = findfirst(i -> resname(ress[i]) in ("PRO", "NPRO", "CPRO") && i > 1, 1:nres)
        @test prolineidx !== nothing
        @test !bp.phirotatable[prolineidx]

        # Map backbone dihedrals to columns.
        idx = 0
        psicol = Dict{Int,Int}()
        phicol = Dict{Int,Int}()
        for i = 1:nres-1
            psicol[i] = (idx += 1)
            bp.phirotatable[i+1] && (phicol[i+1] = (idx += 1))
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

    if testad
        @testset "Differentiability" begin
            function makef(chain)  # for copy/paste command line inferrability
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
