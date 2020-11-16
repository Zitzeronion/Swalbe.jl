using Swalbe
using Documenter

makedocs(;
    modules=[Swalbe],
    authors="Zitzeronion <physiknerd@gmail.com>",
    repo="https://github.com/Ziteronion/Swalbe.jl/blob/{commit}{path}#L{line}",
    sitename="Swalbe.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Ziteronion.github.io/Swalbe.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Ziteronion/Swalbe.jl",
)
