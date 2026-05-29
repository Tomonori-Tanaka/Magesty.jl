# Sunny.jl round-trip validation for the `sce_to_sunny` exporter.
#
# Sunny is a heavy dependency, so this lives in a dedicated environment
# (test/sunny/Project.toml) and runs only via `make test-sunny`; it is NOT part
# of `make test-all`. The check evaluates the generated script's system-building
# section in an actual Sunny session and confirms `energy(sys)` reproduces
# `predict_energy(model, sd) - model.j0` for random spin configurations — the
# correctness gate for the SCE → Sunny conversion (signs, normalization, the
# s = 1 single-ion factor, bond offsets, and primitive unfolding).

using Test
using Magesty
using Sunny
using Random
using LinearAlgebra: norm, inv
using StaticArrays: SVector

const _INTEG = joinpath(@__DIR__, "..", "integration")
_fix(parts...) = joinpath(_INTEG, parts...)

# Evaluate only the system-construction portion of a generated script (before the
# magnetic-structure block, which would randomize and minimize the state).
function _eval_system(script::AbstractString)
	head = split(script, "# ---- Magnetic structure")[1]
	m = Module(:GenSunnyScript)
	Core.eval(m, :(using Sunny))
	Core.eval(m, Meta.parseall(head))
	return Core.eval(m, :sys)
end

# Build the Sunny system from a generated script, map the model's supercell atoms
# onto Sunny sites (reshaping to the supercell for the primitive route), and
# return (is_primitive, max energy error over `n` random configs).
function _roundtrip(model::SCEModel; n::Int = 20)
	script = sce_to_sunny(model)
	primitive = occursin("Cell route: primitive", script)
	sys = _eval_system(script)

	L = Matrix{Float64}(model.basis.structure.supercell.lattice_vectors)
	xf = model.basis.structure.supercell.x_frac
	nat = size(xf, 2)

	if primitive
		prim = Magesty._sunny_build_primitive(model).prim
		Lpi = inv(Matrix(prim.latvecs))
		target = reshape_supercell(sys, Matrix(prim.reshape_matrix))
		sites = [position_to_site(target, Lpi * (L * SVector{3, Float64}(xf[:, a]))) for a = 1:nat]
	else
		target = sys
		sites = [(1, 1, 1, a) for a = 1:nat]
	end
	@test length(unique(sites)) == nat   # site mapping is a bijection

	rng = MersenneTwister(20260529)
	maxerr = 0.0
	for _ = 1:n
		sd = randn(rng, 3, nat)
		for c = 1:nat
			sd[:, c] ./= norm(sd[:, c])
		end
		for a = 1:nat
			set_dipole!(target, sd[:, a], sites[a])
		end
		maxerr = max(maxerr, abs(energy(target) - (Magesty.predict_energy(model, sd) - model.j0)))
	end
	return primitive, maxerr
end

@testset "Sunny round-trip" begin
	# Clean models (cutoff < L/2) take the unfolded primitive route; chain is the
	# multi-cell case (cutoff = L/2, ntran = 2). Tolerance allows for the script's
	# 12-significant-digit coupling literals.
	@testset "primitive route: $(f[1])" for f in
											(("dimer", "dimer.xml"),
		("dimer", "dimer_dmi.xml"),
		("chain", "chain.xml"))
		model = Magesty.load(SCEModel, _fix(f...))
		primitive, err = _roundtrip(model)
		@test primitive
		@test err < 1e-9
	end

	# The minimal 2×2×2 fixtures have interactions beyond half the supercell and
	# fall back to the exact (folded) explicit route; fept also exercises the
	# single-ion path.
	@testset "explicit route: $(f[1])" for f in
										   (("febcc_2x2x2_pm", "scecoeffs.xml"),
		("fept_tetragonal_2x2x2", "scecoeffs.xml"))
		model = Magesty.load(SCEModel, _fix(f...))
		primitive, err = _roundtrip(model)
		@test !primitive
		@test err < 1e-8
	end
end
