#!/usr/bin/env julia
#
# Load EMBSET, vary the number of training structures n along start:step:end, and
# write tables (and optional plots) for whether fitted jphi and energy prediction
# error converge toward the all-data reference fit.
#
# From the Magesty.jl repo root:
#   julia --project=. tools/check_convergence_embset.jl -i input.toml -s system.jld2 -e EMBSET.dat --n-range 5:5:100
#
# Plots.jl is required for PNG output; without it, only TSV files are written.
# Also writes {prefix}_jphi_wide.tsv and {prefix}_jphi_long.tsv (per-coefficient vs. n_train).

using ArgParse
using JLD2
using LinearAlgebra
using Printf
using Random
using Statistics
using TOML

include(joinpath(@__DIR__, "..", "src", "Magesty.jl"))
using .Magesty

function parse_n_range(s::AbstractString)
	parts = split(s, ':')
	if length(parts) == 2
		a, b = parse.(Int, strip.(parts))
		return a:1:b
	elseif length(parts) == 3
		a, st, b = parse.(Int, strip.(parts))
		return a:st:b
	else
		throw(ArgumentError("--n-range must be \"start:end\" or \"start:step:end\": $s"))
	end
end

function energy_rmse_full(
	X_energy::Matrix{Float64},
	j0::Float64,
	jphi::Vector{Float64},
	energies::Vector{Float64},
)::Float64
	pred = X_energy[:, 2:end] * jphi .+ j0
	sqrt(mean((pred .- energies) .^ 2))
end

function _fmt_tsv_e(x::Float64)
	@sprintf("%.15e", x)
end

# Called via Base.invokelatest after `eval(Main, :(using Plots))` so `Main.Plots` is visible
# (same-function `Plots.plot` after eval hits a world-age UndefVarError).
function _save_convergence_embset_figure(
	ns_p::AbstractVector,
	je_p::AbstractVector,
	er_p::AbstractVector,
	png_path::AbstractString,
)
	Plt = Base.getproperty(Main, :Plots)
	p1 = Plt.plot(
		ns_p,
		je_p;
		marker = :circle,
		xlabel = "Training set size n",
		ylabel = "||jphi(n) - jphi(full)|| / ||jphi(full)||",
		title = "jphi convergence (relative L2 vs. full-data reference)",
		legend = false,
	)
	p2 = Plt.plot(
		ns_p,
		er_p;
		marker = :square,
		xlabel = "Training set size n",
		ylabel = "RMSE (eV)",
		title = "Energy prediction RMSE over all N structures",
		legend = false,
	)
	Plt.plot(p1, p2; layout = (2, 1), size = (600, 700))
	Plt.savefig(png_path)
	return nothing
end

function main()
	s = ArgParseSettings(
		description = """
		Scan training-set size n and report convergence of jphi (vs. all-data reference)
		and energy RMSE over all structures.""",
	)
	@add_arg_table s begin
		"--input", "-i"
		help = "Path to input.toml"
		arg_type = String
		required = true
		"--system", "-s"
		help = "system.jld2 (@save system). Omit to build System from input.toml"
		arg_type = String
		default = ""
		"--embset", "-e"
		help = "EMBSET file (default: regression.datafile from input.toml)"
		arg_type = String
		default = ""
		"--n-range", "-n"
		help = "Training counts to evaluate: \"start:end\" or \"start:step:end\" (e.g. 10:10:200)"
		arg_type = String
		required = true
		"--out", "-o"
		help = "Output prefix: .tsv, _jphi_wide.tsv, _jphi_long.tsv, and .png"
		arg_type = String
		default = "convergence_embset"
		"--shuffle"
		help = "Shuffle EMBSET order, then train on the first n structures"
		action = :store_true
		"--seed"
		help = "RNG seed for --shuffle"
		arg_type = Int
		default = 42
		"--no-plot"
		help = "Do not save figure"
		action = :store_true
	end

	args = parse_args(ARGS, s)
	input_path = args["input"]
	system_path = args["system"]
	embset_path = strip(args["embset"])
	out_prefix = args["out"]
	nr = parse_n_range(args["n-range"])

	input = TOML.parsefile(input_path)
	input_dir = dirname(abspath(input_path))
	if isempty(embset_path)
		embset_path = input["regression"]["datafile"]::String
	end
	if !isfile(embset_path)
		alt = isabspath(embset_path) ? "" : joinpath(input_dir, embset_path)
		if !isempty(alt) && isfile(alt)
			embset_path = alt
		else
			error("EMBSET not found: $embset_path (also tried relative to input.toml directory)")
		end
	end

	system = if !isempty(strip(system_path))
		p = strip(system_path)
		if !isfile(p)
			error("System file not found: $p")
		end
		d = JLD2.load(p)
		if !haskey(d, "system")
			error("JLD2 file has no key \"system\" (save with @save \"file.jld2\" system)")
		end
		d["system"]::System
	else
		System(input; verbosity = false)
	end

	all_configs = read_embset(embset_path)
	N = length(all_configs)
	if N == 0
		error("EMBSET contains no structures")
	end

	perm = collect(1:N)
	if args["shuffle"]
		Random.seed!(args["seed"])
		shuffle!(perm)
	end
	ordered_configs = all_configs[perm]

	energies = Float64[c.energy for c in ordered_configs]
	observed_torque_full = [c.torques for c in ordered_configs]

	# Reference fit: all N structures (design matrices live in sc_ref.optimize)
	@info "Fitting all $N structures (reference jphi / j0)…"
	sc_ref = SpinCluster(system, input, ordered_configs; verbosity = false)
	opt_ref = sc_ref.optimize
	X_full = opt_ref.design_matrix_energy
	torque_row_block = 3 * system.structure.supercell.num_atoms
	j0_ref, jphi_ref = Magesty.get_j0_jphi(sc_ref)
	norm_jphi_ref = norm(jphi_ref)
	denom_jphi = norm_jphi_ref > 0 ? norm_jphi_ref : 1.0

	rmse_ref_full = energy_rmse_full(X_full, j0_ref, jphi_ref, energies)
	@info @sprintf("Reference: ||jphi||=%.6g, all-data energy RMSE=%.6g eV", norm_jphi_ref, rmse_ref_full)

	n_list = Int[n for n in nr if 1 <= n <= N]
	n_total = length(n_list)
	if n_total == 0
		@warn "No n in --n-range lies in 1:$N; subset loop skipped."
	else
		@info @sprintf(
			"Subset scan: %d step(s), n ∈ [%d, %d] (full set size N=%d)",
			n_total,
			first(n_list),
			last(n_list),
			N,
		)
	end

	ns = Int[]
	jphi_l2 = Float64[]
	jphi_rel = Float64[]
	energy_rmse = Float64[]
	status = String[]
	j0_list = Float64[]
	coef_rows = Vector{Float64}[]

	M = length(jphi_ref)

	cfg_opt = Magesty.Config4Optimize(input)
	estimator = Magesty.ElasticNet(alpha = cfg_opt.alpha, lambda = cfg_opt.lambda)
	fit_w = cfg_opt.weight

	for (k, n) in enumerate(n_list)
		try
			@views Xe_n = opt_ref.design_matrix_energy[1:n, :]
			@views Xt_n = opt_ref.design_matrix_torque[1:(n * torque_row_block), :]
			@views e_n = energies[1:n]
			@views tau_n = observed_torque_full[1:n]
			j0_n, jphi_n = Magesty.Optimize._fit_sce_model_internal(
				Xe_n,
				Xt_n,
				e_n,
				tau_n,
				estimator,
				fit_w,
			)
			dj = jphi_n .- jphi_ref
			l2 = norm(dj)
			rel = l2 / denom_jphi
			rmse = energy_rmse_full(X_full, j0_n, jphi_n, energies)
			push!(ns, n)
			push!(jphi_l2, l2)
			push!(jphi_rel, rel)
			push!(energy_rmse, rmse)
			push!(status, "ok")
			push!(j0_list, j0_n)
			push!(coef_rows, copy(jphi_n))
			@info @sprintf(
				"[%d/%d] n=%d: ||jphi||=%.6g, all-data energy RMSE=%.6g eV (rel. ||Δjphi|| vs full=%.6g)",
				k,
				n_total,
				n,
				norm(jphi_n),
				rmse,
				rel,
			)
		catch err
			@warn "Fit failed for n=$n" exception = (err, catch_backtrace())
			push!(ns, n)
			push!(jphi_l2, NaN)
			push!(jphi_rel, NaN)
			push!(energy_rmse, NaN)
			push!(status, "fail")
			push!(j0_list, NaN)
			push!(coef_rows, fill(NaN, M))
			@info @sprintf("[%d/%d] n=%d: fit failed (see warning above)", k, n_total, n)
		end
	end

	tsv_path = out_prefix * ".tsv"
	open(tsv_path, "w") do io
		println(
			io,
			join(
				[
					"n",
					"jphi_l2_error_vs_full",
					"jphi_relative_l2_error",
					"energy_rmse_all_configs_eV",
					"status",
				],
				'\t',
			),
		)
		for i in eachindex(ns)
			println(
				io,
				join(
					[
						string(ns[i]),
						@sprintf("%.12e", jphi_l2[i]),
						@sprintf("%.12e", jphi_rel[i]),
						@sprintf("%.12e", energy_rmse[i]),
						status[i],
					],
					'\t',
				),
			)
		end
	end
	@info "Wrote table: $tsv_path"

	jphi_wide_path = out_prefix * "_jphi_wide.tsv"
	open(jphi_wide_path, "w") do io
		println(
			io,
			"# n_train=-1: reference fit on all N structures (N=$N). Other rows: fit on first n_train configs.",
		)
		jphi_hdr = ["jphi_salc_$i" for i in 1:M]
		println(io, join(vcat(["n_train", "status", "j0"], jphi_hdr), '\t'))
		ref_line = String[string(-1), "reference", _fmt_tsv_e(j0_ref)]
		for j in 1:M
			push!(ref_line, _fmt_tsv_e(jphi_ref[j]))
		end
		println(io, join(ref_line, '\t'))
		for i in eachindex(ns)
			line = String[string(ns[i]), status[i], _fmt_tsv_e(j0_list[i])]
			row = coef_rows[i]
			for j in 1:M
				push!(line, _fmt_tsv_e(row[j]))
			end
			println(io, join(line, '\t'))
		end
	end
	@info "Wrote per-coefficient table (wide): $jphi_wide_path"

	jphi_long_path = out_prefix * "_jphi_long.tsv"
	open(jphi_long_path, "w") do io
		println(
			io,
			join(
				["n_train", "salc_index", "jphi", "jphi_ref_full", "delta_vs_ref"],
				'\t',
			),
		)
		for salc in 1:M
			refv = jphi_ref[salc]
			println(
				io,
				join(
					[-1, salc, _fmt_tsv_e(refv), _fmt_tsv_e(refv), _fmt_tsv_e(0.0)],
					'\t',
				),
			)
		end
		for i in eachindex(ns)
			row = coef_rows[i]
			for salc in 1:M
				v = row[salc]
				refv = jphi_ref[salc]
				dlt = isnan(v) ? NaN : v - refv
				println(
					io,
					join(
						[
							ns[i],
							salc,
							_fmt_tsv_e(v),
							_fmt_tsv_e(refv),
							_fmt_tsv_e(dlt),
						],
						'\t',
					),
				)
			end
		end
	end
	@info "Wrote per-coefficient table (long): $jphi_long_path"

	if !args["no-plot"]
		plot_ok = false
		if Base.find_package("Plots") === nothing
			@warn "Plots not found; install with `using Pkg; Pkg.add(\"Plots\")` to save PNG."
		else
			try
				Base.eval(Main, :(using Plots))
				mask = .!isnan.(energy_rmse)
				ns_p = ns[mask]
				je_p = jphi_rel[mask]
				er_p = energy_rmse[mask]
				png_path = out_prefix * ".png"
				Base.invokelatest(
					_save_convergence_embset_figure,
					ns_p,
					je_p,
					er_p,
					png_path,
				)
				@info "Saved figure: $png_path"
				plot_ok = true
			catch err
				@warn "Saving figure with Plots failed" exception = (err, catch_backtrace())
			end
		end
		if !plot_ok && Base.find_package("Plots") !== nothing
			@info "Figure not saved; see table at $tsv_path."
		end
	end

	return nothing
end

main()
