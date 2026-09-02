using DihedralParametrization
using Documenter

DocMeta.setdocmeta!(DihedralParametrization, :DocTestSetup, :(using DihedralParametrization); recursive=true)

makedocs(;
    modules=[DihedralParametrization],
    authors="Tim Holy <tim.holy@gmail.com> and contributors",
    sitename="DihedralParametrization.jl",
    checkdocs=:public,
    format=Documenter.HTML(;
        canonical="https://JuliaStructBio.github.io/DihedralParametrization.jl",
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
    repo="github.com/JuliaStructBio/DihedralParametrization.jl",
    devbranch="main",
)
