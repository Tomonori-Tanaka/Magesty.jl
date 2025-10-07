using Documenter
using Magesty

DocMeta.setdocmeta!(Magesty, :DocTestSetup, :(using Magesty); recursive=true)

makedocs(
    sitename = "Magesty.jl",
    modules  = [Magesty],
    remotes  = nothing,
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://Tomonori-Tanaka.github.io/Magesty.jl",
        assets = String[],
    ),
    pages    = [
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Tutorial" => "tutorial.md",
        "Examples" => "examples.md",
        "Tools" => "tools.md",
    ],
    warnonly = true,
    checkdocs = :exports,
    doctest = true,
)