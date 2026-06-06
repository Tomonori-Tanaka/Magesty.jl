"""
	IncarIO

Reader and writer for VASP `INCAR` text files. Parameters are kept in an
`OrderedDict{Symbol, Any}` so file order is preserved across a parse/write
round trip.

This is the VASP-specific I/O adapter used by the spin-configuration sampler;
it is distinct from `VaspIO`, which parses the XML `vasprun.xml`.
"""
module IncarIO

using DataStructures: OrderedDict
using Printf

export parse_incar, write_incar

"""
	parse_incar(filename::AbstractString)::OrderedDict{Symbol, Any}

Parse a VASP `INCAR` file into an `OrderedDict{Symbol, Any}`, preserving the
order in which parameters appear.

# Arguments
- `filename::AbstractString`: path to the INCAR file.

# Returns
- `OrderedDict{Symbol, Any}`: parameter name => value, in file order.

# Examples
```julia
params = parse_incar("INCAR")
params[:MAGMOM]   # Vector{Float64}
```

# Notes
- Empty lines and comment lines (starting with `#` or `!`) are ignored.
- Multiple parameters on one line may be separated by semicolons.
- Backslash (`\\`) continues a value onto the next line.
- A bare `KEY = "` opens a multi-line (heredoc) string value, closed by a line
  containing only `"`.
- `RWIGS`, `MAGMOM`, and `M_CONSTR` are always returned as `Vector{Float64}`;
  other space-separated numeric lists become vectors, single numbers become
  scalars (`Int` when integral, else `Float64`).
- Logical values (`.TRUE.`/`T`/`.FALSE.`/`F`) become `Bool`.
- `n*v` repeat patterns are expanded (e.g. `3*0` becomes `0 0 0`).
"""
function parse_incar(filename::AbstractString)::OrderedDict{Symbol, Any}
	if !isfile(filename)
		throw(ArgumentError("INCAR file not found at: $filename"))
	end

	params = OrderedDict{Symbol, Any}()
	current_line = ""
	in_heredoc = false
	heredoc_key = ""
	heredoc_value = ""

	open(filename, "r") do io
		for raw_line in eachline(io)
			line = strip(raw_line)
			if !in_heredoc
				isempty(line) && continue
				startswith(line, '#') && continue
				startswith(line, '!') && continue
			end

			# Heredoc handling.
			if !in_heredoc && occursin(r"^\w+\s*=\s*\"$", line)
				in_heredoc = true
				key_match = match(r"^(\w+)\s*=", line)
				key_match === nothing && continue  # guaranteed by occursin above
				heredoc_key = Symbol(something(key_match.captures[1]))
				heredoc_value = ""
				continue
			elseif in_heredoc && line == "\""
				in_heredoc = false
				params[heredoc_key] = strip(heredoc_value)
				continue
			elseif in_heredoc
				heredoc_value *= line * "\n"
				continue
			end

			# Line continuation with a trailing backslash.
			line_rstrip = rstrip(line)
			if endswith(line_rstrip, "\\")
				current_line *= line_rstrip[1:end-1] * " "
				continue
			else
				current_line *= line * " "
			end

			for param in split(current_line, ';')
				param = strip(param)
				isempty(param) && continue

				parts = split(param, '='; limit = 2)
				length(parts) == 2 || continue
				param_name = Symbol(strip(parts[1]))

				# Strip trailing comments from the value.
				value_part = strip(split(parts[2], '#'; limit = 2)[1])
				value_part = strip(split(value_part, '!'; limit = 2)[1])
				value_part = strip(replace(value_part, r"[\r\n]+" => " "))

				# Drop a trailing unit (e.g. "500 eV").
				value = if occursin(r"\d+\s*[a-zA-Z]+$", value_part)
					strip(split(value_part, r"\s+[a-zA-Z]+$")[1])
				else
					value_part
				end

				# Space-separated numeric lists.
				if occursin(r"^[\d\s\.\-\*]+$", value)
					expanded = _expand_repeats(value)
					numbers = filter(!isempty, split(expanded, r"\s+"))
					if all(occursin(r"^-?\d+\.?\d*$", n) for n in numbers)
						if param_name in (:RWIGS, :MAGMOM, :M_CONSTR)
							params[param_name] = [parse(Float64, n) for n in numbers]
						elseif length(numbers) == 1
							params[param_name] = _parse_number(numbers[1])
						else
							params[param_name] = [_parse_number(n) for n in numbers]
						end
						continue
					end
				end

				# Scalars, logicals, and bare strings.
				if occursin(r"^-?\d+\.?\d*$"i, value)
					params[param_name] = _parse_number(value)
				elseif occursin(r"^\.?true\.?$"i, value) || uppercase(value) == "T"
					params[param_name] = true
				elseif occursin(r"^\.?false\.?$"i, value) || uppercase(value) == "F"
					params[param_name] = false
				else
					params[param_name] = value
				end
			end

			current_line = ""
		end
	end

	if haskey(params, :MAGMOM) && haskey(params, :M_CONSTR)
		magmom_length = length(params[:MAGMOM])
		m_constr_length = length(params[:M_CONSTR])
		if magmom_length != m_constr_length
			@warn "MAGMOM and M_CONSTR have different lengths" magmom_length m_constr_length
		end
	end

	return params
end

"""
	write_incar(filename::AbstractString, params::OrderedDict{Symbol, Any};
	            wrap_vectors::Bool = false)::Nothing

Write `params` to a VASP `INCAR` file, preserving insertion order.

`MAGMOM` / `M_CONSTR` are formatted as 3-component vectors separated by a double
space. With `wrap_vectors = true`, `MAGMOM` / `M_CONSTR` wrap with a trailing
backslash every 8 vectors and other vectors wrap every 30 elements.

# Arguments
- `filename::AbstractString`: output path.
- `params::OrderedDict{Symbol, Any}`: parameters to write.

# Keyword arguments
- `wrap_vectors::Bool = false`: wrap long vectors onto continuation lines.

# Returns
- `Nothing`.

# Examples
```julia
write_incar("INCAR.out", params; wrap_vectors = true)
```
"""
function write_incar(
	filename::AbstractString,
	params::OrderedDict{Symbol, Any};
	wrap_vectors::Bool = false,
)::Nothing
	open(filename, "w") do io
		for (param_name, value) in params
			if value isa AbstractVector
				_write_vector(io, param_name, value, wrap_vectors)
			elseif value isa Bool
				println(io, "$param_name = $(value ? "TRUE" : "FALSE")")
			else
				println(io, "$param_name = $value")
			end
		end
	end
	return nothing
end

# Format a single MAGMOM/M_CONSTR atom block "x y z" with 9 decimals.
function _format_spin_vector(value::AbstractVector, i::Integer)::String
	j = 3 * (i - 1)
	return @sprintf("%.9f %.9f %.9f", value[j+1], value[j+2], value[j+3])
end

function _write_vector(
	io::IO,
	param_name::Symbol,
	value::AbstractVector,
	wrap_vectors::Bool,
)
	is_spin = param_name in (:MAGMOM, :M_CONSTR)
	if wrap_vectors && is_spin
		num_atoms = length(value) ÷ 3
		print(io, "$param_name = ")
		for i = 1:num_atoms
			if i > 1 && (i - 1) % 8 == 0
				print(io, "  \\\n  ")
			end
			vec_str = _format_spin_vector(value, i)
			print(io, i < num_atoms ? "$vec_str  " : vec_str)
		end
		println(io)
	elseif wrap_vectors
		first_chunk = true
		for chunk in Iterators.partition(value, 30)
			if first_chunk
				print(io, "$param_name = $(join(chunk, " "))")
			else
				print(io, "  $(join(chunk, " "))")
			end
			println(io, length(chunk) == 30 ? " \\" : "")
			first_chunk = false
		end
	elseif is_spin
		num_atoms = length(value) ÷ 3
		vec_strings = [_format_spin_vector(value, i) for i = 1:num_atoms]
		println(io, "$param_name = $(join(vec_strings, "  "))")
	else
		println(io, "$param_name = $(join(value, " "))")
	end
	return nothing
end

# Parse "12" as Int, "1.5"/"1." as Float64.
function _parse_number(s::AbstractString)::Union{Int, Float64}
	return occursin(r"^-?\d+$", s) ? parse(Int, s) : parse(Float64, s)
end

# Expand VASP repeat syntax "n*v" into n space-separated copies of v, e.g.
# "3*0.0" -> "0.0 0.0 0.0". Tokens without a "*" are left unchanged.
function _expand_repeats(s::AbstractString)::String
	pat = r"(\d+)\*(-?\d+\.?\d*)"
	return replace(s, pat => function (token)
		m = match(pat, token)
		m === nothing && return String(token)  # only matched tokens are passed
		count = parse(Int, something(m.captures[1]))
		value = String(something(m.captures[2]))
		return join(fill(value, count), " ")
	end)
end

end # module IncarIO
