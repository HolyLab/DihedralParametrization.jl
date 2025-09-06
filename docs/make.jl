using ProteinOptimization
using Documenter

DocMeta.setdocmeta!(ProteinOptimization, :DocTestSetup, :(using ProteinOptimization); recursive=true)

makedocs(;
    modules=[ProteinOptimization],
    authors="Tim Holy <tim.holy@gmail.com> and contributors",
    sitename="ProteinOptimization.jl",
    format=Documenter.HTML(;
        canonical="https://HolyLab.github.io/ProteinOptimization.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/HolyLab/ProteinOptimization.jl",
    devbranch="main",
)
