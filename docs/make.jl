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
		assets = ["assets/codeblock_filenames.css"],
		mathengine = Documenter.MathJax3(),
		edit_link = "main",
		repolink = "https://github.com/Tomonori-Tanaka/Magesty.jl",
	),
	pages = [
		"Home" => "index.md",
		"Installation" => "installation.md",
		"API Reference" => "api.md",
		"Internal API" => "api_internal.md",
		"Tutorial" => "tutorial.md",
		"Input Keys" => "input_keys.md",
		"Examples" => "examples.md",
		"Tools" => "tools.md",
		hide(
			"Theoretical Background" => "theory/index.md",
			[
				"theory/overview.md",
				"theory/spherical_harmonics.md",
				"theory/angular_momentum_coupling.md",
				"theory/symmetry_adaptation.md",
				"theory/orbit_clusters.md",
				"theory/folded_tensor.md",
				"theory/design_matrix_and_fitting.md",
				"technical_notes.md",
			],
		),
		hide(
			"Tips" => "tips/index.md",
			[
				"tips/penalty_parameter_dependence.md",
				"tips/rwigs_dependence.md",
				"tips/mfa_analysis.md",
			],
		),
	],
	warnonly = true,
	checkdocs = :exports,
	doctest = true,
)

deploydocs(
	repo = "github.com/Tomonori-Tanaka/Magesty.jl",
	devbranch = "main",
)
