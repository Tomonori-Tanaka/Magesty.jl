using Test
using Magesty
using Magesty.Fitting: build_design_matrix_energy, build_design_matrix_torque
using Magesty.SpinConfigs: read_embset
using DelimitedFiles

# Regression harness for the design-matrix algorithmic restructuring.
#
# Each restructuring stage relaxes the floating-point reduction order
# differently. Reference matrices are captured into
# `test/regression_fixtures/` and committed to the repo. This testset
# compares the current `build_design_matrix_*` output against those
# references at the configured tolerance.
#
# The tolerance defaults below match the post-restructuring floating-
# point budget documented in the cluster-major torque theory page at
# `docs/src/theory/cluster_major_torque.md`. Callers needing a tighter
# bound (e.g. before a refactor expected to be bit-identical) can
# override via `MAGESTY_REGRESSION_RTOL_*` env vars.

const FIXTURES_DIR = joinpath(@__DIR__, "..", "regression_fixtures")

# Each entry: (label, input_toml, embset, n_configs_cap). Must mirror the
# capture script.
const FIXTURES = [
	("fege_2x2x2",
		joinpath(@__DIR__, "..", "integration", "fege_2x2x2", "input.toml"),
		joinpath(@__DIR__, "..", "integration", "fege_2x2x2", "EMBSET"),
		5),
	("fept_tetragonal_2x2x2",
		joinpath(@__DIR__, "..", "integration", "fept_tetragonal_2x2x2", "input.toml"),
		joinpath(@__DIR__, "..", "integration", "fept_tetragonal_2x2x2", "EMBSET"),
		5),
]

# Default tolerances correspond to the post-D row of the per-stage table in
# the spec. The matrices captured at M0 baseline are bit-identical to the
# current implementation, so these defaults will pass trivially today and
# tighten naturally as each stage lands its own commit.
const DEFAULT_RTOL_ENERGY = 1e-12
const DEFAULT_ATOL_ENERGY = 1e-13
const DEFAULT_RTOL_TORQUE = 1e-10
const DEFAULT_ATOL_TORQUE = 1e-12

function _env_float(name::String, fallback::Float64)::Float64
	v = get(ENV, name, nothing)
	v === nothing && return fallback
	return parse(Float64, v)
end

function _read_reference(path::String)::Matrix{Float64}
	# Strip comment lines (the capture script writes 4 header lines starting
	# with '#') before handing off to readdlm.
	lines = String[]
	open(path, "r") do io
		for line in eachline(io)
			startswith(line, "#") && continue
			push!(lines, line)
		end
	end
	io = IOBuffer(join(lines, "\n"))
	return readdlm(io, '\t', Float64)
end

@testset "Design-matrix equivalence (spec 260524)" begin
	rtol_e = _env_float("MAGESTY_REGRESSION_RTOL_ENERGY", DEFAULT_RTOL_ENERGY)
	atol_e = _env_float("MAGESTY_REGRESSION_ATOL_ENERGY", DEFAULT_ATOL_ENERGY)
	rtol_t = _env_float("MAGESTY_REGRESSION_RTOL_TORQUE", DEFAULT_RTOL_TORQUE)
	atol_t = _env_float("MAGESTY_REGRESSION_ATOL_TORQUE", DEFAULT_ATOL_TORQUE)

	for (label, input_path, embset_path, ncfg_cap) in FIXTURES
		ref_e_path = joinpath(FIXTURES_DIR, label * "_energy.tsv")
		ref_t_path = joinpath(FIXTURES_DIR, label * "_torque.tsv")
		if !(isfile(ref_e_path) && isfile(ref_t_path))
			@warn "Reference fixture missing; skipping" label ref_e_path ref_t_path
			continue
		end

		@testset "$label" begin
			basis = SCEBasis(input_path; verbosity = false)
			spinconfigs_all = read_embset(embset_path)
			ncfg = min(ncfg_cap, length(spinconfigs_all))
			spinconfigs = spinconfigs_all[1:ncfg]

			salc_list = basis.salcbasis.salc_list
			num_atoms = basis.structure.supercell.num_atoms
			symmetry = basis.symmetry

			X_E = build_design_matrix_energy(salc_list, spinconfigs, symmetry; verbosity = false)
			X_T = build_design_matrix_torque(salc_list, spinconfigs, num_atoms, symmetry; verbosity = false)

			ref_E = _read_reference(ref_e_path)
			ref_T = _read_reference(ref_t_path)

			@test size(X_E) == size(ref_E)
			@test size(X_T) == size(ref_T)
			@test isapprox(X_E, ref_E; rtol = rtol_e, atol = atol_e)
			@test isapprox(X_T, ref_T; rtol = rtol_t, atol = atol_t)
		end
	end
end
