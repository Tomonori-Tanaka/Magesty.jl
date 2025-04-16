#!/usr/bin/env julia
# Sample program to convert a VASP structure file to a TOML file

using ArgParse
using TOML
using OrderedCollections
using Printf

"""
    parse_vasp(filename::String) -> NamedTuple

Read and parse the VASP structure file.
Assumes the file follows a common POSCAR format:
  Line 1: Comment
  Line 2: Scaling factor
  Lines 3-5: Lattice vectors
  Line 6: Either element symbols or the number of atoms in each species
  (If line 6 contains symbols, then line 7 will contain the number of atoms)
  Next line: coordinate type (Direct or Cartesian)
  Following lines: atomic positions (first 3 entries are x, y, z coordinates)
"""
function parse_vasp(filename::String)
    if !isfile(filename)
        throw(ErrorException("File not found: $filename"))
    end

    open(filename, "r") do io
        lines = readlines(io)
        
        # Check if file has a minimum number of lines for valid content
        if length(lines) < 9
            throw(ErrorException("Input file appears to be too short to be a valid VASP structure file"))
        end

        # Extract header information
        comment = strip(lines[1])
        scaling = try
            parse(Float64, strip(lines[2]))
        catch e
            throw(ErrorException("Invalid scaling factor in line 2: $(lines[2])"))
        end
        
        # Parse lattice vectors
        lattice_vectors = Vector{Vector{Float64}}(undef, 3)
        for i in 1:3
            try
                lattice_vectors[i] = [parse(Float64, token) for token in split(strip(lines[i+2]))]
                if length(lattice_vectors[i]) != 3
                    throw(ErrorException("Invalid lattice vector in line $(i+2): $(lines[i+2])"))
                end
            catch e
                throw(ErrorException("Error parsing lattice vector in line $(i+2): $(lines[i+2])"))
            end
        end
        
        # Process species and number information
        element_symbols = String[]
        numbers = Int[]
        coord_line_index = 0
        tokens = split(strip(lines[6]))
        
        try
            # Try parsing line 6 as integers (number of atoms)
            parsed_ints = [parse(Int, token) for token in tokens]
            
            if all(x -> x > 0, parsed_ints)
                # Line 6 contains number of atoms, but element symbols are not provided
                numbers = parsed_ints
                element_symbols = fill("X", length(numbers))
                coord_line_index = 7
            else
                throw(ErrorException("Invalid number of atoms in line 6: $(lines[6])"))
            end
        catch
            # Line 6 contains element symbols; then line 7 must contain the number of atoms
            element_symbols = tokens
            try
                numbers = [parse(Int, token) for token in split(strip(lines[7]))]
                if !all(x -> x > 0, numbers)
                    throw(ErrorException("Invalid number of atoms in line 7: $(lines[7])"))
                end
                coord_line_index = 8
            catch e
                throw(ErrorException("Error parsing number of atoms in line 7: $(lines[7])"))
            end
        end

        # Parse coordinate mode
        coordinate_mode = lowercase(strip(lines[coord_line_index]))
        if coordinate_mode ∉ ["direct", "cartesian"]
            throw(ErrorException("Invalid coordinate mode: $coordinate_mode"))
        end

        # Read atomic positions
        n_atoms = sum(numbers)
        start_pos = coord_line_index + 1
        positions = Vector{Vector{Float64}}(undef, n_atoms)
        
        for i in 1:n_atoms
            line_idx = start_pos + i - 1
            if line_idx > length(lines)
                throw(ErrorException("Unexpected end of file while reading atomic positions"))
            end
            
            try
                tokens = split(strip(lines[line_idx]))
                if length(tokens) < 3
                    throw(ErrorException("Invalid position format in line $line_idx"))
                end
                positions[i] = [parse(Float64, token) for token in tokens[1:3]]
            catch e
                throw(ErrorException("Error parsing position in line $line_idx: $(lines[line_idx])"))
            end
        end
        
        # Return parsed data as a NamedTuple
        return (
            comment = comment,
            scaling = scaling,
            lattice_vectors = lattice_vectors,
            coordinate_mode = coordinate_mode,
            element_symbols = element_symbols,
            numbers = numbers,
            positions = positions
        )
    end
end

"""
    create_toml(data::NamedTuple, output_file::String) -> Nothing

Creates a TOML file from parsed VASP structure data.

# Arguments
- `data::NamedTuple`: Parsed VASP structure data
- `output_file::String`: Path to the output TOML file

# Throws
- `ErrorException` if the file cannot be written
"""
function create_toml(data::NamedTuple, output_file::String)
    try
        # Prepare the TOML data structure
        toml_data = OrderedDict{String, Any}()
        
        # [general] section
        general = OrderedDict{String, Any}()
        general["name"] = data.comment
        general["nat"] = sum(data.numbers)
        general["kd"] = data.element_symbols
        general["periodicity"] = [true, true, true]
        general["j_zero_thr"] = 1e-10
        toml_data["general"] = general
        
        # [symmetry] section
        symmetry = OrderedDict{String, Any}()
        symmetry["tolerance"] = 1e-5
        toml_data["symmetry"] = symmetry
        
        # [interaction] section
        interaction = OrderedDict{String, Any}()
        interaction["nbody"] = 2
        
        # interaction.lmax
        lmax = OrderedDict{String, Any}()
        for element in data.element_symbols
            lmax[element] = [0, 1]
        end
        interaction["lmax"] = lmax
        
        # interaction.cutoff
        cutoff = OrderedDict{String, Any}()
        for i in 1:length(data.element_symbols)
            for j in i:length(data.element_symbols)
                pair = join(sort([data.element_symbols[i], data.element_symbols[j]]), "-")
                cutoff[pair] = [-1, -1]
            end
        end
        interaction["cutoff"] = cutoff
        toml_data["interaction"] = interaction
        
        # [regression] section
        regression = OrderedDict{String, Any}()
        regression["weight"] = 0.0
        regression["datafile"] = "EMBSET"
        toml_data["regression"] = regression
        
        # [structure] section
        structure = OrderedDict{String, Any}()
        structure["scale"] = data.scaling
        structure["lattice"] = data.lattice_vectors
        
        # Create kd_list and position arrays
        kd_list = Int[]
        positions = Vector{Vector{Float64}}()
        
        # Assign species numbers (1 for first species, 2 for second, etc.)
        atom_index = 1
        for (i, count) in enumerate(data.numbers)
            for j in 1:count
                push!(kd_list, i)
                push!(positions, data.positions[atom_index])
                atom_index += 1
            end
        end
        
        structure["kd_list"] = kd_list
        structure["position"] = positions
        toml_data["structure"] = structure
        
        # Write the TOML data to the output file with formatting
        open(output_file, "w") do io
            # Write each section with custom formatting
            write(io, "[general]\n")
            for (key, value) in general
                if key == "kd"
                    write(io, "kd = [\"$(join(value, "\", \""))\"]\n")
                elseif key == "periodicity"
                    write(io, "periodicity = [$(join(value, ", "))]\n")
                elseif key == "name"
                    write(io, "name = \"$value\"\n")
                else
                    write(io, "$key = $value\n")
                end
            end
            write(io, "\n")
            
            write(io, "[symmetry]\n")
            for (key, value) in symmetry
                write(io, "$key = $value\n")
            end
            write(io, "\n")
            
            write(io, "[interaction]\n")
            write(io, "nbody = $(interaction["nbody"])\n\n")
            
            write(io, "[interaction.lmax]\n")
            for (element, value) in lmax
                write(io, "$element = [$(join(value, ", "))]\n")
            end
            write(io, "\n")
            
            write(io, "[interaction.cutoff]\n")
            for (pair, value) in cutoff
                write(io, "$pair = [$(join(value, ", "))]\n")
            end
            write(io, "\n")
            
            write(io, "[regression]\n")
            for (key, value) in regression
                if key == "datafile"
                    write(io, "$key = \"$value\"\n")
                else
                    write(io, "$key = $value\n")
                end
            end
            write(io, "\n")
            
            write(io, "[structure]\n")
            write(io, "scale = $(structure["scale"])\n")
            write(io, "lattice = [\n")
            for vec in structure["lattice"]
                write(io, "        [ $(join(vec, ", ")) ],\n")
            end
            write(io, "]\n")
            
            write(io, "kd_list = [\n")
            # 20原子ごと、または数値が変わったときに改行
            current_line = Int[]
            current_value = kd_list[1]
            
            for (i, value) in enumerate(kd_list)
                if value != current_value || length(current_line) >= 20
                    if !isempty(current_line)
                        write(io, "        $(join(current_line, ", ")),\n")
                    end
                    current_line = Int[]
                    current_value = value
                end
                push!(current_line, value)
            end
            
            # 最後の行を出力
            if !isempty(current_line)
                write(io, "        $(join(current_line, ", ")),\n")
            end
            write(io, "]\n")
            
            # Format position array with aligned decimal points
            write(io, "position = [\n")
            for pos in positions
                # Format each coordinate with fixed width
                formatted_pos = map(x -> rpad(@sprintf("%.15f", x), 20), pos)
                write(io, "        [ $(join(formatted_pos, ", ")) ],\n")
            end
            write(io, "]\n")
        end
        println("TOML file generated: $output_file")
    catch e
        throw(ErrorException("Error writing TOML file: $(e.msg)"))
    end
end

"""
    main() -> Nothing

Parses command-line arguments and executes the conversion.

# Throws
- `ErrorException` if command-line arguments are invalid
"""
function main()
    # Set up argument parser
    s = ArgParseSettings(
        description = "Convert VASP structure file to TOML format",
        version = "1.0.0",
        add_version = true
    )
    
    @add_arg_table s begin
        "--input", "-i"
        help = "Input VASP structure file (e.g., POSCAR)"
        arg_type = String
        required = true

        "--output", "-o"
        help = "Output TOML file"
        arg_type = String
        default = "input.toml"
    end
    
    try
        args = parse_args(s)
        input_file = args["input"]
        output_file = args["output"]
        
        # Parse the VASP file
        parsed_data = parse_vasp(input_file)
        
        # Generate the TOML file
        create_toml(parsed_data, output_file)
    catch e
        if isa(e, ArgParseError)
            println(stderr, e.msg)
            exit(1)
        else
            rethrow(e)
        end
    end
end

# Execute the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end