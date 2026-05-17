#!/usr/bin/env julia
# Symmetry benchmark.
# Times `Symmetries.Symmetry(structure, options)` — the spglib call plus
# the `map_sym` / `map_sym_inv` permutation tables — across the three
# committed integration fixtures, so the cost of the symmetry pre-pass
# is visible as a separate line in `.claude/bench_log.md`.
#
# Usage:
#   julia --project=bench bench/benchmark_symmetry.jl

using Magesty
using Magesty.Structures: Structure
using Magesty.Symmetries: Symmetry
using Magesty.InputSpecs: parse_toml_inputs
using Statistics
using Printf
using TOML

const FIXTURES = [
	"febcc_2x2x2_pm",
	"fept_tetragonal_2x2x2",
	"fege_2x2x2",
]
const NTRIALS = 5

# Same shape as bench_b1_design_matrix.bench_one: one warm-up call, then
# `NTRIALS` `@timed` samples for wall-time / bytes, plus a second pass
# with `@allocations` because `@timed` does not report the allocation
# count directly.
function bench_one(label::String, f::Function)
	f()  # warm-up
	times_s = Float64[]
	bytes_b = Float64[]
	for _ in 1:NTRIALS
		stats = @timed f()
		push!(times_s, stats.time)
		push!(bytes_b, stats.bytes)
	end
	counts = Float64[]
	for _ in 1:NTRIALS
		push!(counts, Float64(@allocations f()))
	end
	@printf("%-44s  min=%.4fs  med=%.4fs  bytes(med)=%.2e  allocs(med)=%d\n",
		label,
		minimum(times_s),
		median(times_s),
		median(bytes_b),
		Int(median(counts)))
end

function load_specs(fixture::AbstractString)
	toml_path = joinpath(
		@__DIR__, "..", "test", "integration", fixture, "input.toml",
	)
	input = TOML.parsefile(toml_path)
	system, _interaction, options = parse_toml_inputs(input)
	return system, options
end

function main()
	println("=== Symmetry benchmark (5 trials each) ===")
	for fixture in FIXTURES
		system, options = load_specs(fixture)
		structure = Structure(system; verbosity = false)
		println()
		println("Fixture: ", fixture,
			" (num_atoms = ", structure.supercell.num_atoms, ")")
		bench_one("Symmetry(structure, options)",
			() -> Symmetry(structure, options; verbosity = false))
	end
end

main()
