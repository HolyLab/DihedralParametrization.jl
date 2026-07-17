using DihedralParametrization
using BioStructures
using StaticArrays
using LinearAlgebra
using ForwardDiff
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

    @testset "Bond lengths and angles are invariant under every dihedral" begin
        # The parametrization's defining contract: reconstructing at any
        # dihedral vector must leave every bond length and every bond angle
        # fixed at its encoded value, since only the dihedrals themselves are
        # free parameters. This must hold across residue boundaries too --
        # C(i)-N(i+1), O(i)-C(i)-N(i+1), C(i)-N(i+1)-CA(i+1), and
        # C(i)-N(i+1)-H(i+1) -- since a same-residue build reference for H(i)
        # or O(i) would let it co-rotate with φ(i) or ψ(i) instead of staying
        # fixed. A test that only looked at same-residue geometry would pass
        # vacuously.
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

        # Connectivity comes from BioStructures' residue templates, which give
        # each residue's intra-residue bonds plus the atoms carrying bonds to
        # adjacent residues. Deriving it from interatomic distances instead
        # would need element-dependent cutoffs (S-H is ~1.34 Å against C-H's
        # ~1.09 Å) and would still have to special-case geminal H pairs.
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

        # Confirm the topology actually reaches the cross-residue quartets
        # where the documented violations live.
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
        # Proline's ring nitrogen carries no amide H, so C(i)-N(i+1)-H(i+1)
        # exists only where residue i+1 has one.
        @test nCNH == count(i -> haskey(ress[i].atoms, "H"), 2:nres)

        # ∂(bond length)/∂θ_k and ∂(bond angle)/∂θ_k for every bond, angle,
        # and dihedral index k, from a single ForwardDiff pass over the
        # reconstruction. Any nonzero entry is a violation of the invariance
        # contract. These derivatives are a detector, not a magnitude: at a
        # planar amide a misplaced H or O sits at a local extremum of its
        # spurious rotation, so the angle derivative can read ~0.04 deg/rad
        # while the atom itself moves ~1 Å/rad. Any nonzero entry is a bug;
        # its size says little about how far the atom strays.
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
        # Proline's ring (N-CA-CB-CG-CD-N) has rigid bond lengths and angles
        # and so has no internal mobility: chi1, chi2, and phi are all fixed
        # by the ring rather than free parameters. This is the only source of
        # residue-count-dependent DOF loss relative to a naive per-residue count.
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
