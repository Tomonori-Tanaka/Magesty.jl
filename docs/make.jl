import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.develop(Pkg.PackageSpec(path = joinpath(@__DIR__, "..")))
using Magesty
using Documenter

DocMeta.setdocmeta!(Magesty, :DocTestSetup, :(using Magesty); recursive = true)

makedocs(
	sitename = "Magesty.jl",
	modules = [Magesty],
	remotes = nothing,
	format = Documenter.HTML(
		prettyurls = get(ENV, "CI", "false") == "true",
		canonical = "https://Tomonori-Tanaka.github.io/Magesty.jl",
		assets = String[],
		mathengine = Documenter.MathJax3(),
		edit_link = "main",
		repolink = "https://github.com/Tomonori-Tanaka/Magesty.jl",
	),
	pages = [
		"Home" => "index.md",
		"API Reference" => "api.md",
		"Tutorial" => "tutorial.md",
		"Examples" => "examples.md",
		"Tools" => "tools.md",
		"Technical Notes" => "technical_notes.md",
		"Tips" => [
			"Overview" => "tips/index.md",
			"Penalty term dependence" => "tips/penalty_parameter_dependence.md",
			"RWIGS dependence" => "tips/rwigs_dependence.md",
		],
	],
	warnonly = true,
	checkdocs = :none,
	doctest = true,
)

deploydocs(
	repo = "github.com/Tomonori-Tanaka/Magesty.jl",
	devbranch = "main",
)
