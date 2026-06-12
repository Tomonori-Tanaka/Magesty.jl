import Pkg
Pkg.activate(@__DIR__)
Pkg.develop(Pkg.PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.instantiate()
using Magesty
using Documenter
using Literate

DocMeta.setdocmeta!(Magesty, :DocTestSetup, :(using Magesty); recursive = true)

# Generate the worked-example pages (markdown for Documenter, plus a
# downloadable script and notebook) from their Literate.jl sources.
let literate_dir = joinpath(@__DIR__, "literate"),
	examples_out = joinpath(@__DIR__, "src", "examples")

	for src in ("case1_bcc_fe.jl",)
		path = joinpath(literate_dir, src)
		# The hand-written front half (prerequisites, file preparation) is plain
		# markdown; only the runnable workflow comes from the Literate source.
		# Prepend the intro to the generated page, but keep the downloadable
		# script and notebook limited to the workflow code.
		intro_path = joinpath(literate_dir, replace(src, ".jl" => "_intro.md"))
		intro = isfile(intro_path) ? read(intro_path, String) : ""
		Literate.markdown(
			path, examples_out;
			documenter = true,
			postprocess = md -> intro * "\n\n" * md,
		)
		Literate.script(path, examples_out)
		Literate.notebook(path, examples_out; execute = false)
	end
end

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
		"API Examples" => "api_examples.md",
		"Tutorial" => "tutorial.md",
		"Input Keys" => "input_keys.md",
		"Examples" => [
			"examples/index.md",
			"examples/case1_bcc_fe.md",
		],
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
				"theory/sh_cache.md",
				"theory/cluster_major_torque.md",
				"theory/design_matrix_and_fitting.md",
				"technical_notes.md",
			],
		),
		hide(
			"Tips" => "tips/index.md",
			[
				"tips/penalty_parameter_dependence.md",
				"tips/rwigs_dependence.md",
				"tips/mfa_sampling.md",
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
