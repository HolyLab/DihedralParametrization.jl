# DihedralParametrization

<!--
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HolyLab.github.io/DihedralParametrization.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HolyLab.github.io/DihedralParametrization.jl/dev/)
-->
[![Build Status](https://github.com/HolyLab/DihedralParametrization.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/HolyLab/DihedralParametrization.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/HolyLab/DihedralParametrization.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/HolyLab/DihedralParametrization.jl)

DihedralParametrization represents [protein structures](https://en.wikipedia.org/wiki/Protein_structure)
using their rotatable-bond [dihedral angles](https://en.wikipedia.org/wiki/Dihedral_angle#In_stereochemistry).
These include the backbone angles `ϕ` and `ψ` and side-chain `χ` angles.

If you are only interested in the backbone degrees of freedom, see [Backboner](https://github.com/MurrellGroup/Backboner.jl) instead.


## Installation

DihedralParametrization is not registered. Install it with

```
pkg> add https://github.com/HolyLab/DihedralParametrization.jl
```

## Pre-requisites

Input must be a complete structure loaded from a CIF or PDB file by
[BioStructures](https://github.com/BioJulia/BioStructures.jl). It must satisfy
these constraints:

- Hydrogens use the Amber naming convention. [ChimeraX](https://www.cgl.ucsf.edu/chimerax/)
  can add them.
- Histidine is disambiguated as
  [HID/HIE/HIP](https://ambermd.org/Questions/HIS.html), and the amino- and
  carboxyl-terminal residues should have "N" and "C" prepended to their residue
  names, respectively. BioStructure's `specializeresnames!(chain)` should
  handle these names after hydrogens have been added.

## Demo

The example uses an AlphaFold2 structure from `test/data` with hydrogens added.

```julia
julia> using DihedralParametrization, BioStructures

julia> cd(joinpath(pkgdir(DihedralParametrization), "test", "data"))

julia> struc = read("AF-M3YHX5-F1-model_v4_hydrogens.cif", MMCIFFormat)
MolecularStructure AF-M3YHX5-F1-model_v4_hydrogens.cif with 1 models, 1 chains (A), 112 residues, 1833 atoms

julia> specializeresnames!(struc)
MolecularStructure AF-M3YHX5-F1-model_v4_hydrogens.cif with 1 models, 1 chains (A), 112 residues, 1833 atoms

julia> A = struc["A"]
Chain A with 112 residues, 0 other molecules, 1833 atoms

julia> A[1]    # check that we've renamed the amino-terminal residue
Residue 1:A with name NMET, 19 atoms

julia> bp, dihedrals = bondparametrization(A);

julia> bp
BondParametrization with 1833 atoms and 112 residues

julia> summary(dihedrals)
"562-element Vector{Float64}"
```

`bp` stores bond lengths, bond angles, and fixed dihedral angles. `dihedrals`
stores the rotatable dihedral angles in radians.

Pass a different vector with entries from `-π` to `π` to generate another
configuration:

```julia
julia> randdh = 2π * (rand(length(dihedrals)) .- 0.5);

julia> X = atomcoordinates(bp, randdh, A);

julia> summary(X)
"1833-element Vector{SVector{3, Float64}}"

julia> B = buildchain(A, bp, X)
Chain A with 112 residues, 0 other molecules, 1833 atoms
```

Plot the original and generated structures:

```julia
using GLMakie, ProtPlot
fig = Figure(size=(1000, 500))
ax1 = LScene(fig[1, 1]; show_axis=false)
ax2 = LScene(fig[1, 2]; show_axis=false)
ribbon!(ax1, A)
ribbon!(ax2, B)
```

![structure comparison](docs/src/assets/randchain.png)
