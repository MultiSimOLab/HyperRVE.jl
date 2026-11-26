using HyperRVE
using Documenter

DocMeta.setdocmeta!(HyperRVE, :DocTestSetup, :(using HyperRVE); recursive=true)

makedocs(;
    modules=[HyperRVE],
    authors="MultiSimO Lab",
    sitename="HyperRVE.jl",
    format=Documenter.HTML(;
        canonical="https://multisimolab.github.io/HyperRVE.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md", "User guide" => "guide.md", "Examples" => "examples.md",
           "Migration from v1" => "migration.md"],
    warnonly=[:missing_docs],
)

deploydocs(; repo="github.com/MultiSimOLab/HyperRVE.jl", devbranch="main")
