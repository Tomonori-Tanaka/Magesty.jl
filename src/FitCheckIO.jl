# Plain-text writers for fit-quality inspection.
#
# `write_energies` and `write_torques` dump observed (DFT) versus predicted
# (SCE) values to whitespace-separated text files. The files are consumed by
# the `FitCheck_energy.py` / `FitCheck_torque.py` visualization scripts under
# `tools/`. The writers are a thin I/O layer: predictions go through the
# `predict_energy` / `predict_torque` verbs unchanged, and observed values are
# read off the `SpinConfig`s.

# --- Argument resolution ------------------------------------------------

# Predictions are produced by a SCEModel; a SCEFit is reduced to one.
_resolve_predictor(model::SCEModel)::SCEModel = model
_resolve_predictor(f::SCEFit)::SCEModel = SCEModel(f)

# Observed values are read off SpinConfigs. The three accepted data forms
# all reduce to a vector of SpinConfig. The result is only iterated, so a
# vector passed directly is returned without copying.
_resolve_spinconfigs(data::SCEDataset) = data.spinconfigs
_resolve_spinconfigs(data::AbstractVector{SpinConfig}) = data
_resolve_spinconfigs(embset_path::AbstractString) = read_embset(embset_path)

# Per-atom element symbols of the predictor's structure, in atom order.
function _element_names(model::SCEModel)::Vector{String}
	structure = model.basis.structure
	return String[
		structure.kd_name[kd] for kd in structure.supercell.kd_int_list
	]
end

# --- File writers -------------------------------------------------------

# Write the energy file: a two-line header followed by one row per
# configuration with columns `index  DFT_Energy  SCE_Energy` (eV).
function _write_energy_file(
	filename::AbstractString,
	observed::AbstractVector{<:Real},
	predicted::AbstractVector{<:Real},
)::Nothing
	n_configs = length(observed)
	length(predicted) == n_configs || throw(ArgumentError(
		"observed ($n_configs) and predicted ($(length(predicted))) " *
		"energy counts differ"))
	idx_width = max(ndigits(n_configs), 1)
	row_fmt = Printf.Format(" %$(idx_width)d    % 15.10e    % 15.10e\n")
	try
		open(filename, "w") do io
			println(io, "# data index,    DFT_Energy,    SCE_Energy")
			println(io, "# unit of energy is eV")
			for i = 1:n_configs
				Printf.format(io, row_fmt, i, observed[i], predicted[i])
			end
		end
	catch e
		@error "Failed to write energies to $filename" exception =
			(e, catch_backtrace())
		rethrow(e)
	end
	return nothing
end

# Write the torque file: a two-line header followed by a `# data index: N`
# block per configuration, each block holding one row per atom with columns
# `atom_index  element  DFT_xyz  SCE_xyz` (eV).
function _write_torque_file(
	filename::AbstractString,
	observed::AbstractVector{<:AbstractMatrix{<:Real}},
	predicted::AbstractVector{<:AbstractMatrix{<:Real}},
	element_names::AbstractVector{<:AbstractString},
)::Nothing
	n_configs = length(observed)
	length(predicted) == n_configs || throw(ArgumentError(
		"observed ($n_configs) and predicted ($(length(predicted))) " *
		"torque counts differ"))
	num_atoms = length(element_names)
	idx_width = max(ndigits(num_atoms), 1)
	element_width = isempty(element_names) ? 1 : maximum(length, element_names)
	row_fmt = Printf.Format(
		" %$(idx_width)d %$(element_width)s  % 15.10e   % 15.10e   " *
		"% 15.10e    % 15.10e   % 15.10e   % 15.10e\n")
	try
		open(filename, "w") do io
			println(io,
				"# atom index,    element,   DFT_torque_x,    " *
				"DFT_torque_y,    DFT_torque_z,    SCE_torque_x,    " *
				"SCE_torque_y,    SCE_torque_z")
			println(io, "# unit of torque is eV")
			for c = 1:n_configs
				obs = observed[c]
				pred = predicted[c]
				size(obs) == (3, num_atoms) || throw(ArgumentError(
					"observed torque matrix for configuration $c has size " *
					"$(size(obs)), expected (3, $num_atoms)"))
				size(pred) == (3, num_atoms) || throw(ArgumentError(
					"predicted torque matrix for configuration $c has size " *
					"$(size(pred)), expected (3, $num_atoms)"))
				println(io, "# data index: $c")
				for a = 1:num_atoms
					Printf.format(io, row_fmt,
						a, element_names[a],
						obs[1, a], obs[2, a], obs[3, a],
						pred[1, a], pred[2, a], pred[3, a])
				end
			end
		end
	catch e
		@error "Failed to write torques to $filename" exception =
			(e, catch_backtrace())
		rethrow(e)
	end
	return nothing
end

# --- Public API ---------------------------------------------------------

"""
	write_energies(f::SCEFit, filename::AbstractString = "energy_list.txt")
	write_energies(predictor::Union{SCEModel, SCEFit},
	               data::Union{SCEDataset, AbstractVector{SpinConfig}, AbstractString},
	               filename::AbstractString)

Write observed (DFT) versus predicted (SCE) energies to a text file for
fit-quality inspection. The output is consumed by the `FitCheck_energy.py`
visualization script under `tools/`.

The two-argument form evaluates `f` against its own training dataset; the
lone string argument is always the output path. The three-argument form
evaluates any predictor against an explicit dataset, e.g. a held-out
validation or test set, and requires an explicit output path so the data
and output arguments cannot be confused.

# Arguments
- `f::SCEFit`: A fitted model; evaluated against `f.dataset`.
- `predictor::Union{SCEModel, SCEFit}`: The predictor producing SCE energies.
- `data::Union{SCEDataset, AbstractVector{SpinConfig}, AbstractString}`: The
  configurations supplying observed energies. An `AbstractString` is read as
  an EMBSET file path.
- `filename::AbstractString`: Output path. Defaults to `"energy_list.txt"`
  only in the two-argument `SCEFit` form.

# Returns
- `nothing`. The file at `filename` is created (overwritten if it exists).

The file has a two-line comment header followed by one row per
configuration: `data_index  DFT_Energy  SCE_Energy`. Energies are in the
unit of the training data (typically eV).

# Examples
```julia
f = fit(SCEFit, dataset, Ridge(lambda = 1e-4))
write_energies(f)                                       # -> energy_list.txt
write_energies(f, "train_E.txt")                        # training set
write_energies(f, test_dataset, "test_E.txt")           # held-out set
write_energies(SCEModel(f), "EMBSET", "E.txt")          # model + EMBSET path
```
"""
write_energies(f::SCEFit, filename::AbstractString = "energy_list.txt")::Nothing =
	write_energies(SCEModel(f), f.dataset, filename)

function write_energies(
	predictor::Union{SCEModel,SCEFit},
	data::Union{SCEDataset,AbstractVector{SpinConfig},AbstractString},
	filename::AbstractString,
)::Nothing
	model = _resolve_predictor(predictor)
	spinconfigs = _resolve_spinconfigs(data)
	observed = Float64[sc.energy for sc in spinconfigs]
	predicted = predict_energy(model, spinconfigs)
	return _write_energy_file(filename, observed, predicted)
end

"""
	write_torques(f::SCEFit, filename::AbstractString = "torque_list.txt")
	write_torques(predictor::Union{SCEModel, SCEFit},
	              data::Union{SCEDataset, AbstractVector{SpinConfig}, AbstractString},
	              filename::AbstractString)

Write observed (DFT) versus predicted (SCE) per-atom torques to a text file
for fit-quality inspection. The output is consumed by the `FitCheck_torque.py`
visualization script under `tools/`.

The two-argument form evaluates `f` against its own training dataset; the
lone string argument is always the output path. The three-argument form
evaluates any predictor against an explicit dataset, e.g. a held-out
validation or test set, and requires an explicit output path so the data
and output arguments cannot be confused.

# Arguments
- `f::SCEFit`: A fitted model; evaluated against `f.dataset`.
- `predictor::Union{SCEModel, SCEFit}`: The predictor producing SCE torques.
- `data::Union{SCEDataset, AbstractVector{SpinConfig}, AbstractString}`: The
  configurations supplying observed torques. An `AbstractString` is read as
  an EMBSET file path.
- `filename::AbstractString`: Output path. Defaults to `"torque_list.txt"`
  only in the two-argument `SCEFit` form.

# Returns
- `nothing`. The file at `filename` is created (overwritten if it exists).

The file has a two-line comment header followed by a `# data index: N` block
per configuration. Each block holds one row per atom:
`atom_index  element  DFT_torque_{x,y,z}  SCE_torque_{x,y,z}`. Torque
components are in the unit of the training data (typically eV).

# Examples
```julia
f = fit(SCEFit, dataset, Ridge(lambda = 1e-4))
write_torques(f)                                # -> torque_list.txt
write_torques(f, "train_T.txt")                 # training set
write_torques(f, test_dataset, "test_T.txt")    # held-out set
write_torques(SCEModel(f), "EMBSET", "T.txt")   # model + EMBSET path
```
"""
write_torques(f::SCEFit, filename::AbstractString = "torque_list.txt")::Nothing =
	write_torques(SCEModel(f), f.dataset, filename)

function write_torques(
	predictor::Union{SCEModel,SCEFit},
	data::Union{SCEDataset,AbstractVector{SpinConfig},AbstractString},
	filename::AbstractString,
)::Nothing
	model = _resolve_predictor(predictor)
	spinconfigs = _resolve_spinconfigs(data)
	observed = Matrix{Float64}[sc.torques for sc in spinconfigs]
	predicted = predict_torque(model, spinconfigs)
	element_names = _element_names(model)
	return _write_torque_file(filename, observed, predicted, element_names)
end

# --- GCV diagnostic writers ---------------------------------------------

"""
	write_gcv_lambda(path::GCVLambdaPath,
	                 filename::AbstractString = "gcv_lambda.txt")

Write a ridge GCV penalty sweep ([`gcv_lambda`](@ref)) to a text file for
inspection or plotting. The output is consumed by the `FitCheck_gcv_lambda.py`
script under `tools/`.

The file has a comment header followed by one row per penalty:
`lambda  gcv  dof`. The GCV score is in the weighted-objective unit (not eV²),
and `dof` is the effective degrees of freedom `tr(H)`. The selected
`lambda_best` and the torque weight are recorded in the header.

# Arguments
- `path::GCVLambdaPath`: The sweep result.
- `filename::AbstractString`: Output path. Defaults to `"gcv_lambda.txt"`.

# Returns
- `nothing`. The file at `filename` is created (overwritten if it exists). Any
  filesystem error is logged and re-thrown.

# Examples
```julia
path = gcv_lambda(dataset, 10.0 .^ (-6:0.5:0))
write_gcv_lambda(path, "gcv_lambda.txt")
```
"""
function write_gcv_lambda(
	path::GCVLambdaPath,
	filename::AbstractString = "gcv_lambda.txt",
)::Nothing
	n = length(path.lambdas)
	row_fmt = Printf.Format(" % 15.10e    % 15.10e    % 15.10e\n")
	try
		open(filename, "w") do io
			println(io, "# ridge GCV penalty sweep (gcv_lambda)")
			println(io, "# torque_weight = ", path.torque_weight,
				",  lambda_best = ", path.lambda_best)
			println(io, "# GCV is in the weighted-objective unit (not eV^2)")
			println(io, "# lambda,    GCV,    effective_dof")
			for i = 1:n
				Printf.format(io, row_fmt,
					path.lambdas[i], path.gcv_scores[i], path.dof[i])
			end
		end
	catch e
		@error "Failed to write GCV lambda path to $filename" exception =
			(e, catch_backtrace())
		rethrow(e)
	end
	return nothing
end

"""
	write_gcv_learning_curve(curve::GCVSizeCurve,
	                         filename::AbstractString = "gcv_learning_curve.txt")

Write a data-sufficiency GCV learning curve ([`gcv_learning_curve`](@ref)) to a
text file for inspection or plotting. The output is consumed by the
`FitCheck_gcv_learning_curve.py` script under `tools/`.

The file has a comment header followed by one row per training-set size:
`size  gcv_mean  gcv_std`, where the mean and standard deviation are taken over
the random subset draws at that size. The GCV score is in the weighted-objective
unit (not eV²). The estimator, torque weight, repeats, and seed are recorded in
the header.

# Arguments
- `curve::GCVSizeCurve`: The sweep result.
- `filename::AbstractString`: Output path. Defaults to
  `"gcv_learning_curve.txt"`.

# Returns
- `nothing`. The file at `filename` is created (overwritten if it exists). Any
  filesystem error is logged and re-thrown.

# Examples
```julia
curve = gcv_learning_curve(dataset, Ridge(lambda = 1e-4); repeats = 8)
write_gcv_learning_curve(curve, "gcv_learning_curve.txt")
```
"""
function write_gcv_learning_curve(
	curve::GCVSizeCurve,
	filename::AbstractString = "gcv_learning_curve.txt",
)::Nothing
	n = length(curve.sizes)
	size_width = max(maximum(ndigits, curve.sizes; init = 1), 4)
	row_fmt = Printf.Format(" %$(size_width)d    % 15.10e    % 15.10e\n")
	try
		open(filename, "w") do io
			println(io, "# data-sufficiency GCV learning curve (gcv_learning_curve)")
			println(io, "# estimator = ", curve.estimator)
			println(io, "# torque_weight = ", curve.torque_weight,
				",  repeats = ", curve.repeats, ",  seed = ", curve.seed)
			println(io, "# GCV is in the weighted-objective unit (not eV^2)")
			println(io, "# size,    GCV_mean,    GCV_std")
			for i = 1:n
				Printf.format(io, row_fmt,
					curve.sizes[i], curve.gcv_mean[i], curve.gcv_std[i])
			end
		end
	catch e
		@error "Failed to write GCV learning curve to $filename" exception =
			(e, catch_backtrace())
		rethrow(e)
	end
	return nothing
end
