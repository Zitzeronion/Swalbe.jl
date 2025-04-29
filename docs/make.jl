using Documenter, Swalbe

# doctest(Swalbe; manual = false)
DocMeta.setdocmeta!(Swalbe, :DocTestSetup, :(using Swalbe); recursive=true)
makedocs(
    modules=[Swalbe],
    authors="Zitzeronion <physiknerd@gmail.com>",
    repo="https://github.com/Zitzeronion/Swalbe.jl",
    sitename="Swalbe.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://Ziteronion.github.io/Swalbe.jl",
        assets=String[],
        size_threshold = nothing
    ),
    # checkdocs = :exports,
    pages=[
        "Home" => "index.md",
        "Tutorials" => "tutorials.md",
        "Functions" => "functions.md",
    ],
)

deploydocs(
    repo = "github.com/Zitzeronion/Swalbe.jl.git",
)
