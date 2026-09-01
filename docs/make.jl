using DihedralParametrization
using Documenter

DocMeta.setdocmeta!(DihedralParametrization, :DocTestSetup, :(using DihedralParametrization); recursive=true)

makedocs(;
    modules=[DihedralParametrization],
    authors="Tim Holy <tim.holy@gmail.com> and contributors",
    sitename="DihedralParametrization.jl",
    checkdocs=:all,
    format=Documenter.HTML(;
        canonical="https://HolyLab.github.io/DihedralParametrization.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Derivatives" => "derivatives.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/HolyLab/DihedralParametrization.jl",
    devbranch="main",
)
