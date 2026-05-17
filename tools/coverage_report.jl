#!/usr/bin/env julia
# Per-file test coverage report for src/.
#
# Reads the `.cov` files produced by `julia --code-coverage=user` (set by
# `make test-coverage`), aggregates them with Coverage.jl, and prints a
# table sorted by ascending coverage percentage so the files most in
# need of additional tests appear at the top. Also writes an LCOV file
# at the repository root for downstream consumers (Codecov, VS Code
# Coverage Gutters, etc.).
#
# Usage:
#   julia --project=coverage tools/coverage_report.jl

using Coverage
using Printf

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const SRC_DIR   = joinpath(REPO_ROOT, "src")
const LCOV_PATH = joinpath(REPO_ROOT, "coverage.lcov")

function summarize(fc::FileCoverage)
	covered, coverable = get_summary(fc)
	pct = coverable == 0 ? 100.0 : 100.0 * covered / coverable
	return (
		file = relpath(fc.filename, REPO_ROOT),
		covered = covered,
		coverable = coverable,
		pct = pct,
	)
end

function main()
	coverage = process_folder(SRC_DIR)
	if isempty(coverage)
		println(stderr, "No coverage data found under $(SRC_DIR).")
		println(stderr, "Run `make test-coverage` (or the equivalent) first.")
		exit(1)
	end

	rows = [summarize(fc) for fc in coverage]
	# Files with zero coverable lines (e.g. pure module re-exports) carry
	# no signal; keep them but push to the bottom by treating their pct
	# as 100% (already done in `summarize`).
	sort!(rows; by = r -> (r.pct, -r.coverable, r.file))

	println("Coverage report (ascending by coverage %)")
	println(repeat("=", 70))
	for r in rows
		@printf("%-50s %6.1f%%  (%4d / %4d)\n",
			r.file, r.pct, r.covered, r.coverable)
	end
	println(repeat("-", 70))

	total_covered   = sum(r -> r.covered, rows)
	total_coverable = sum(r -> r.coverable, rows)
	total_pct = total_coverable == 0 ? 100.0 :
		100.0 * total_covered / total_coverable
	@printf("%-50s %6.1f%%  (%4d / %4d)\n",
		"TOTAL", total_pct, total_covered, total_coverable)

	LCOV.writefile(LCOV_PATH, coverage)
	println()
	println("LCOV written to ", relpath(LCOV_PATH, REPO_ROOT))
end

main()
