using RVEs
using Documenter

DocMeta.setdocmeta!(RVEs, :DocTestSetup, :(using RVEs); recursive=true)

makedocs(;
    modules=[RVEs],
    authors="MultiSimO Lab",
    sitename="RVEs.jl",
    format=Documenter.HTML(;
        canonical="https://multisimolab.github.io/RVEs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md", "Migration from v1" => "migration.md"],
    warnonly=[:missing_docs],
)

deploydocs(; repo="github.com/MultiSimOLab/RVEs.jl", devbranch="main")
