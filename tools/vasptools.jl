"""
    VaspTools

A module for handling VASP-related file operations.
"""
module VaspTools

using OrderedCollections
using Printf

export parse_incar, write_incar

"""
    parse_incar(filename::String)::OrderedDict{Symbol, Any}

Parse a VASP INCAR file and return its contents as an OrderedDict.

# Arguments
- `filename::String`: Path to the INCAR file

# Returns
- `OrderedDict{Symbol, Any}`: An OrderedDict containing all key-value pairs from the INCAR file, preserving the order of parameters

# Examples
```julia
julia> params = parse_incar("INCAR")
OrderedDict{Symbol, Any} with 19 entries:
  :SYSTEM  => "AIMD"
  :EDIFF   => 0.0001
  :ISYM    => 0
  :ENCUT   => 500.0
  :LREAL   => "A"
  :PREC    => "N"
  :ISMEAR  => 0
  :SIGMA   => 0.05
  :IBRION  => 0
  :NSW     => 10000
  :LCHARG  => false
  :LWAVE   => false
  :SMASS   => 3.0
  :NBLOCK  => 1
  :KBLOCK  => 50
  :POTIM   => 1.0
  :APACO   => 10
  :NPACO   => 200
  :TEBEG   => 300.0
  :TEEND   => 300.0
  :NPAR    => 4
```

# Notes
- Empty lines and comment lines (starting with # or !) are ignored
- Multiple parameters can be specified in one line using semicolons
- Space-separated lists of numbers are converted to arrays
- Line continuation with backslash (`\\`) is supported
- Multi-line string values (heredoc) are supported using double quotes
- Units (e.g., "eV") are stripped from numeric values
- Numeric values are converted to Float64
- Boolean values are converted to Bool (true/false, T/F, .TRUE./.FALSE.)
- The order of parameters is preserved as they appear in the file
"""
function parse_incar(filename::String)::OrderedDict{Symbol, Any}
	if !isfile(filename)
		throw(ArgumentError("INCAR file not found at: $filename"))
	end

	params = OrderedDict{Symbol, Any}()
	current_line = ""
	in_heredoc = false
	heredoc_key = ""
	heredoc_value = ""
	
	open(filename, "r") do io
		for line in eachline(io)
			# Skip empty lines and comment lines (unless in heredoc)
			line = strip(line)
			if !in_heredoc
				isempty(line) && continue
				startswith(line, '#') && continue
				startswith(line, '!') && continue
			end
			
			# Handle heredoc
			if !in_heredoc && occursin(r"^\w+\s*=\s*\"$", line)
				# Start of heredoc
				in_heredoc = true
				heredoc_key = Symbol(match(r"^(\w+)\s*=", line)[1])
				heredoc_value = ""
				continue
			elseif in_heredoc && line == "\""
				# End of heredoc
				in_heredoc = false
				params[heredoc_key] = strip(heredoc_value)
				continue
			elseif in_heredoc
				# Inside heredoc
				heredoc_value *= line * "\n"
				continue
			end
			
			# Handle line continuation
			if endswith(rstrip(line), '\\')
				current_line *= rstrip(line[1:end-1]) * " "
				continue
			else
				current_line *= line * " "
			end
			
			# Process the complete line (including continued lines)
			for param in split(current_line, ';')
				param = strip(param)
				isempty(param) && continue
				
				# Split key and value, handling comments after the value
				parts = split(param, '=', limit=2)
				if length(parts) == 2
					param_name = Symbol(strip(parts[1]))
					# Split value and comment
					value_part = strip(split(parts[2], '#', limit=2)[1])
					value_part = strip(split(value_part, '!', limit=2)[1])
					
					# Remove units from numeric values
					value = if occursin(r"\d+\s*[a-zA-Z]+$", value_part)
						strip(split(value_part, r"\s+[a-zA-Z]+$")[1])
					else
						value_part
					end
					
					# Handle space-separated lists of numbers
					if occursin(r"^[\d\s\.\-\*]+$", value)
						# Expand patterns like "3*0" to "0 0 0"
						expanded_value = replace(value, r"(\d+)\*([\d\.\-]+)" => s"\1 * \2")
						# Split into individual numbers and parse
						numbers = filter(!isempty, split(expanded_value, r"\s+"))
						if all(occursin(r"^-?\d+\.?\d*$", n) for n in numbers)
							# Always treat RWIGS, MAGMOM, M_CONSTR as vectors
							if param_name in [:RWIGS, :MAGMOM, :M_CONSTR]
								params[param_name] = [parse(Float64, n) for n in numbers]
							else
								# For other keys, use single value if there's only one number
								if length(numbers) == 1
									params[param_name] = parse_number(numbers[1])
								else
									params[param_name] = [parse_number(n) for n in numbers]
								end
							end
							continue
						end
					end
					
					# Convert to appropriate type
					if occursin(r"^-?\d+\.?\d*$"i, value)
						params[param_name] = parse_number(value)
					elseif occursin(r"^\.?true\.?$"i, value) || uppercase(value) == "T"
						params[param_name] = true
					elseif occursin(r"^\.?false\.?$"i, value) || uppercase(value) == "F"
						params[param_name] = false
					else
						params[param_name] = value
					end
				end
			end
			
			# Reset current_line for the next iteration
			current_line = ""
		end
	end

	# Check if both MAGMOM and M_CONSTR are present and have the same length
	if haskey(params, :MAGMOM) && haskey(params, :M_CONSTR)
		magmom_length = length(params[:MAGMOM])
		m_constr_length = length(params[:M_CONSTR])
		if magmom_length != m_constr_length
			@warn "MAGMOM and M_CONSTR have different lengths: MAGMOM has $magmom_length elements, M_CONSTR has $m_constr_length elements"
		end
	end
	
	return params
end

"""
    write_incar(filename::String, params::OrderedDict{Symbol, Any}; wrap_vectors::Bool=false)

Write parameters to a VASP INCAR file.

# Arguments
- `filename::String`: Path to the output INCAR file
- `params::OrderedDict{Symbol, Any}`: Parameters to write to the file
- `wrap_vectors::Bool`: If true, wrap vectors every 30 elements with a backslash

# Examples
```julia
julia> params = OrderedDict(
    :SYSTEM => "AIMD",
    :ENCUT => 500.0,
    :LCHARG => false
)
julia> write_incar("INCAR", params)
julia> write_incar("INCAR", params, wrap_vectors=true)  # Wrap vectors every 30 elements
```
"""
function write_incar(filename::String, params::OrderedDict{Symbol, Any}; wrap_vectors::Bool=false)
	open(filename, "w") do io
		for (param_name, value) in params
			if value isa Vector
				# Handle arrays (e.g., MAGMOM)
				if wrap_vectors
					# Split into chunks of 30 elements
					chunks = Iterators.partition(value, 30)
					first_chunk = true
					for chunk in chunks
						if !first_chunk
							print(io, "  ")  # Add indentation for continuation lines
						end
						if first_chunk
							if param_name in [:MAGMOM, :M_CONSTR]
								print(io, "$param_name = $(join([@sprintf("%.9f", x) for x in chunk], " "))")
							else
								print(io, "$param_name = $(join(chunk, " "))")
							end
						else
							if param_name in [:MAGMOM, :M_CONSTR]
								print(io, "$(join([@sprintf("%.9f", x) for x in chunk], " "))")
							else
								print(io, "$(join(chunk, " "))")
							end
						end
						if length(chunk) == 30
							println(io, " \\")
						else
							println(io)
						end
						first_chunk = false
					end
				else
					if param_name in [:MAGMOM, :M_CONSTR]
						println(io, "$param_name = $(join([@sprintf("%.9f", x) for x in value], " "))")
					else
						println(io, "$param_name = $(join(value, " "))")
					end
				end
			elseif value isa Bool
				# Handle boolean values
				println(io, "$param_name = $(value ? "TRUE" : "FALSE")")
			else
				println(io, "$param_name = $value")
			end
		end
	end
end

"""
    parse_number(s::AbstractString)::Union{Int, Float64}

Parse a string into either an integer or a float, depending on its format.

# Arguments
- `s::AbstractString`: The string to parse

# Returns
- `Union{Int, Float64}`: The parsed number as an integer if possible, otherwise as a float
"""
function parse_number(s::AbstractString)::Union{Int, Float64}
	if occursin(r"^-?\d+$", s)
		return parse(Int, s)
	else
		return parse(Float64, s)
	end
end

end
