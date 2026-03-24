"""
    plot_jphi_cluster_distance.jl

Read `jphi.xml` (Magesty `write_xml` output) and plot each SALC coefficient `jphi`
against the maximum over cluster atom pairs of minimum-image (MIC) distance.

- By default, all SALCs (all N-body terms) are shown.
- Different N (`body`) are drawn as separate series / colors.
- Restrict to selected N with `--bodies 2,3`.

Requires Plots.jl (`using Pkg; Pkg.add("Plots")`). When saving with `-o`, sets
`GKSwstype=100` for GR if unset (headless-friendly).
"""

using LinearAlgebra
using EzXML
using ArgParse
using Printf
using Plots

if !@isdefined(Magesty)
	include(joinpath(@__DIR__, "..", "src", "Magesty.jl"))
end
using .Magesty

"""Minimum-image distance (Å) between atoms i and j; same 27-cell sweep as `plot_jij_atom.jl`."""
function pairwise_mic_distance(st::Magesty.Structure, i::Int, j::Int)::Float64
	i == j && return 0.0
	min_d = Inf
	for cell in 1:27  # center cell + neighbor image cells (Structures.NUM_CELLS)
		st.exist_image[cell] || continue
		d = norm(
			st.x_image_cart[:, i, 1] - st.x_image_cart[:, j, cell],
		)
		min_d = min(min_d, d)
	end
	return min_d
end

"""Maximum of the above over all pairs in `atoms`; single-atom cluster returns 0.0."""
function max_intracluster_mic_distance(st::Magesty.Structure, atoms::Vector{Int})::Float64
	n = length(atoms)
	n < 2 && return 0.0
	m = 0.0
	for ii in 1:(n - 1), jj in (ii + 1):n
		m = max(m, pairwise_mic_distance(st, atoms[ii], atoms[jj]))
	end
	return m
end

function parse_bodies_filter(s::AbstractString)::Union{Nothing, Vector{Int}}
	st = strip(s)
	if isempty(st) || lowercase(st) == "all"
		return nothing
	end
	parts = split(st, ','; keepempty = false)
	return [parse(Int, strip(p)) for p in parts]
end

function load_salc_metadata(xml_path::String)::Dict{Int, Tuple{Int, Vector{Int}}}
	doc = readxml(xml_path)
	system_node = findfirst("//System", doc)
	isnothing(system_node) && throw(ArgumentError("<System> not found in XML: $xml_path"))
	basis_node = findfirst("SCEBasisSet", system_node)
	isnothing(basis_node) && throw(ArgumentError("<SCEBasisSet> not found in XML: $xml_path"))
	out = Dict{Int, Tuple{Int, Vector{Int}}}()
	for salc in findall("SALC", basis_node)
		idx = parse(Int, salc["index"])
		body = parse(Int, salc["body"])
		bnodes = findall("basis", salc)
		isempty(bnodes) && continue
		atoms = parse.(Int, split(bnodes[1]["atoms"]))
		if length(atoms) != body
			@warn "SALC index=$idx: body=$body does not match number of atoms in first basis ($(length(atoms))); using first basis atoms for distance."
		end
		out[idx] = (body, atoms)
	end
	return out
end

function load_jphi_values(xml_path::String)::Vector{Tuple{Int, Float64}}
	doc = readxml(xml_path)
	jnode = findfirst("//JPhi", doc)
	isnothing(jnode) &&
		throw(ArgumentError("<JPhi> not found (use XML written with write_jphi=true): $xml_path"))
	out = Tuple{Int, Float64}[]
	for el in findall("jphi", jnode)
		si = parse(Int, el["salc_index"])
		val = parse(Float64, nodecontent(el))
		push!(out, (si, val))
	end
	sort!(out, by = x -> x[1])
	return out
end

function collect_plot_rows(
	xml_path::String,
	bodies_filter::Union{Nothing, Vector{Int}},
)::Vector{NamedTuple{(:salc_index, :nbody, :r_max, :jphi), Tuple{Int, Int, Float64, Float64}}}
	st = Magesty.Structure(xml_path, verbosity = false)
	meta = load_salc_metadata(xml_path)
	rows = NamedTuple{(:salc_index, :nbody, :r_max, :jphi), Tuple{Int, Int, Float64, Float64}}[]
	for (salc_index, jphi) in load_jphi_values(xml_path)
		if !haskey(meta, salc_index)
			@warn "No <SALC> for salc_index=$salc_index in XML; skipping."
			continue
		end
		body, atoms = meta[salc_index]
		if bodies_filter !== nothing && !(body in bodies_filter)
			continue
		end
		rmax = max_intracluster_mic_distance(st, atoms)
		push!(rows, (; salc_index, nbody = body, r_max = rmax, jphi))
	end
	return rows
end

function plot_rows(
	rows;
	title_str::String = "SCE jφ vs max intracluster distance",
	ylabel_unit::String = "eV",
	output_path::Union{Nothing, String} = nothing,
)
	if isempty(rows)
		@warn "No data to plot."
		return nothing
	end

	dists = [r.r_max for r in rows]
	vals = [r.jphi for r in rows]
	# Group labels: distinct N -> distinct series / colors in Plots
	labels = ["N = $(r.nbody)" for r in rows]

	p = scatter(
		dists,
		vals;
		group = labels,
		markersize = 6,
		markeralpha = 0.85,
		title = title_str,
		xlabel = "Max pairwise MIC distance (Å)",
		ylabel = "jφ ($ylabel_unit)",
		legend = :best,
		grid = true,
		size = (800, 600),
		dpi = 300,
		framestyle = :box,
	)
	hline!(p, [0.0]; color = :black, linestyle = :dash, linewidth = 1, alpha = 0.45, label = "")

	println(@sprintf("%-6s %5s %14s %16s", "salc", "N", "R_max(Å)", "jphi"))
	println("----------------------------------------------")
	for r in sort(rows, by = x -> (x.nbody, x.r_max, x.salc_index))
		println(@sprintf("%-6d %5d %14.6f %16.8e", r.salc_index, r.nbody, r.r_max, r.jphi))
	end
	println("----------------------------------------------")
	println("points: ", length(rows))

	if output_path !== nothing
		savefig(p, output_path)
		println("Saved: ", output_path)
	else
		display(p)
		println("\nPress Enter to exit...")
		readline()
	end
	return p
end

function main()
	s = ArgParseSettings(; description = "Plot jphi vs max intracluster MIC distance from Magesty XML.")
	@add_arg_table s begin
		"input"
		help = "Path to jphi.xml (or other Magesty write_xml output)"
		required = true
		arg_type = String
		"--bodies", "-b"
		help = "Comma-separated N-body (body) values to include; omit or 'all' for all clusters."
		arg_type = String
		default = "all"
		"--output", "-o"
		help = "Output path (PNG, etc.). If set, save and exit without interactive display."
		arg_type = String
		default = ""
		"--title"
		help = "Plot title"
		arg_type = String
		default = "SCE jφ vs max intracluster distance"
	end
	args = parse_args(ARGS, s)
	xml_path = args["input"]
	bodies_f = parse_bodies_filter(args["bodies"])
	outp = String(strip(args["output"]))
	outp = isempty(outp) ? nothing : outp
	# Headless / no display backend: only when saving and GKSwstype not already set
	if outp !== nothing && !haskey(ENV, "GKSwstype")
		ENV["GKSwstype"] = "100"
	end

	rows = collect_plot_rows(xml_path, bodies_f)
	plot_rows(
		rows;
		title_str = args["title"],
		output_path = outp,
	)
	return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
