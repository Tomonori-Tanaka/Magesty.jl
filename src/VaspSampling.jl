# MFA spin-configuration sampling for VASP. Orchestrates the `IncarIO` reader /
# writer and the code-agnostic `MfaSampling` sweep into the `sample_mfa_incar`
# public API.

"""
	sample_mfa_incar(incar_path; variable, start, stop, num_points,
	                 num_samples=1, randomize=false, fix="", uniform_atoms="",
	                 outdir=".", prefix="sample") -> Vector{String}

Sample thermally conditioned spin configurations from a VASP `INCAR` and write
each one to its own `INCAR` file.

The initial spin matrix is read from `MAGMOM` (or `M_CONSTR` if `MAGMOM` is
absent). For every value in an evenly spaced sweep of the control variable,
`num_samples` configurations are drawn with the Mean-Field Approximation
sampler (per-atom directions from a von Mises-Fisher distribution, magnitudes
preserved). Each output file copies all keys from the input `INCAR` and sets
both `MAGMOM` and `M_CONSTR` to the sampled configuration.

# Arguments
- `incar_path::AbstractString`: path to the input `INCAR`.

# Keyword arguments
- `variable::AbstractString` (required): control variable, `"tau"` (scaled
  temperature `T/Tc`, expected in `(0, 1]`) or `"m"` (magnetization, expected
  in `[0, 1)`). Values outside the range are clamped to the corresponding
  ordered/disordered limit.
- `start::Real`, `stop::Real`, `num_points::Integer` (required): the sweep
  values are `range(start, stop; length = num_points)`.
- `num_samples::Integer = 1`: configurations drawn per sweep value.
- `randomize::Bool = false`: apply a random global rotation (quantization-axis
  randomization) to each drawn configuration.
- `fix::AbstractString = ""`: 1-based atom indices kept at their input
  directions (rotated by the same global rotation when `randomize`), e.g.
  `"1-10,12,20-22"`.
- `uniform_atoms::AbstractString = ""`: 1-based atom indices drawn isotropically
  on the sphere instead of from the vMF distribution (same index syntax).
- `outdir::AbstractString = "."`: directory for the output files (created if
  needed).
- `prefix::AbstractString = "sample"`: output file-name prefix. The
  `magesty vasp mfa` command does not expose this; it always uses `"sample"`.

# Returns
- `Vector{String}`: the written file paths, in sweep order. Files are named
  `joinpath(outdir, "<prefix>-NN.INCAR")` where `NN` runs from `1` to
  `num_points * num_samples` in `(point, sample)` order, zero-padded to the
  width of that count.

# Examples
```julia
# 3 temperatures from 0.1 to 0.3, two samples each -> 6 INCAR files.
sample_mfa_incar("INCAR"; variable="tau", start=0.1, stop=0.3,
                 num_points=3, num_samples=2, outdir="samples")
```
"""
function sample_mfa_incar(
	incar_path::AbstractString;
	variable::AbstractString,
	start::Real,
	stop::Real,
	num_points::Integer,
	num_samples::Integer = 1,
	randomize::Bool = false,
	fix::AbstractString = "",
	uniform_atoms::AbstractString = "",
	outdir::AbstractString = ".",
	prefix::AbstractString = "sample",
)::Vector{String}
	variable in ("tau", "m") ||
		throw(ArgumentError("variable must be \"tau\" or \"m\"; got \"$variable\""))

	incar = IncarIO.parse_incar(incar_path)
	magmom = if haskey(incar, :MAGMOM)
		incar[:MAGMOM]
	elseif haskey(incar, :M_CONSTR)
		incar[:M_CONSTR]
	else
		throw(ArgumentError("MAGMOM or M_CONSTR not found in INCAR file: $incar_path"))
	end
	magmom isa AbstractVector ||
		throw(ArgumentError("MAGMOM/M_CONSTR must be a numeric vector; got $(typeof(magmom))"))
	length(magmom) % 3 == 0 ||
		throw(ArgumentError("MAGMOM/M_CONSTR length $(length(magmom)) is not a multiple of 3"))

	n_atoms = length(magmom) ÷ 3
	spin_matrix = reshape(Float64.(magmom), 3, n_atoms)

	fixed_indices = MfaSampling.parse_atom_index_spec(fix; max_index = n_atoms)
	uniform_indices = MfaSampling.parse_atom_index_spec(uniform_atoms; max_index = n_atoms)

	configs = MfaSampling.mfa_sweep(
		spin_matrix;
		variable = variable,
		start = start,
		stop = stop,
		num_points = num_points,
		num_samples = num_samples,
		randomize = randomize,
		fixed_indices = fixed_indices,
		uniform_indices = uniform_indices,
	)

	mkpath(outdir)
	width = ndigits(num_points * num_samples)
	paths = Vector{String}(undef, length(configs))
	for (k, config) in enumerate(configs)
		flat = reshape(config, :)
		out_incar = copy(incar)
		out_incar[:MAGMOM] = flat
		out_incar[:M_CONSTR] = flat
		filename = "$(prefix)-$(lpad(k, width, '0')).INCAR"
		path = joinpath(outdir, filename)
		IncarIO.write_incar(path, out_incar; wrap_vectors = true)
		paths[k] = path
	end
	return paths
end
