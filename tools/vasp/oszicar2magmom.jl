#!/usr/bin/env julia
"""
Extract magnetic moments from the last SCF step in VASP OSZICAR file
and output them in INCAR MAGMOM format (one line).
"""

using ArgParse
using Printf

"""
Parse OSZICAR file and extract magnetic moments for the last SCF step.
"""
function parse_oszicar_last_step(filepath::String, magmom_type::String="MW_int")
    data = Dict{Int, Vector{Float64}}()  # step => [mx1, my1, mz1, mx2, my2, mz2, ...]
    current_step = nothing
    in_magmom_section = false
    pending_data = Vector{Tuple{Int, Float64, Float64, Float64}}()  # (ion_index, mx, my, mz)

    # Pattern to match SCF step markers (3 uppercase letters followed by ':')
    step_pattern = r"^([A-Z]{3}):\s+(\d+)"

    open(filepath, "r") do f
        for line in eachline(f)
            # Detect start of SCF step
            m = match(step_pattern, line)
            if m !== nothing
                # When we encounter a DAV: line, assign pending data to this step
                # The data before DAV: N belongs to step N
                step_number = parse(Int, m.captures[2])
                if !isempty(pending_data)
                    # Sort by ion index and flatten to [mx1, my1, mz1, mx2, my2, mz2, ...]
                    sorted_data = sort(pending_data, by=x -> x[1])
                    magmom_vector = Float64[]
                    for (ion_idx, mx, my, mz) in sorted_data
                        push!(magmom_vector, mx, my, mz)
                    end
                    data[step_number] = magmom_vector
                    empty!(pending_data)
                else
                    # Initialize empty vector if no pending data
                    data[step_number] = Float64[]
                end

                # Set the current step number for next iteration
                current_step = step_number
                in_magmom_section = false
                continue
            end

            # Detect start of magnetic moment section
            if occursin("ion", line) && occursin("MW_int", line) && occursin("M_int", line)
                in_magmom_section = true
                continue
            end

            # Read magnetic moment data
            if in_magmom_section
                parts = split(line)
                if length(parts) >= 7
                    try
                        ion_index = parse(Int, parts[1])

                        if magmom_type == "MW_int"
                            # MW_int x, y, z components (indices 2, 3, 4)
                            mx = parse(Float64, parts[2])
                            my = parse(Float64, parts[3])
                            mz = parse(Float64, parts[4])
                        elseif magmom_type == "M_int"
                            # M_int x, y, z components (indices 5, 6, 7)
                            mx = parse(Float64, parts[5])
                            my = parse(Float64, parts[6])
                            mz = parse(Float64, parts[7])
                        else
                            throw(ArgumentError("Unknown magmom_type: $magmom_type"))
                        end

                        # Store in pending_data
                        push!(pending_data, (ion_index, mx, my, mz))
                    catch e
                        # Skip if numerical conversion fails
                        continue
                    end
                end
            end

            # Detect end of section or next SCF step
            if in_magmom_section && (match(r"^[A-Z]{3}:", line) !== nothing || isempty(strip(line)))
                in_magmom_section = false
            end
        end
    end

    # Assign any remaining pending data to the last step
    if !isempty(pending_data) && current_step !== nothing
        sorted_data = sort(pending_data, by=x -> x[1])
        magmom_vector = Float64[]
        for (ion_idx, mx, my, mz) in sorted_data
            push!(magmom_vector, mx, my, mz)
        end
        data[current_step] = magmom_vector
    end

    return data
end

"""
Format magnetic moment vector as INCAR MAGMOM line.
"""
function format_magmom_line(magmom_vector::Vector{Float64})
    if isempty(magmom_vector)
        return "MAGMOM ="
    end
    
    # Format each value with appropriate precision
    formatted_values = [@sprintf("%.6f", val) for val in magmom_vector]
    return "MAGMOM = " * join(formatted_values, "  ")
end

function main()
    s = ArgParseSettings(
        description="Extract magnetic moments from the last SCF step in OSZICAR file and output in INCAR MAGMOM format",
        epilog="""
Examples:
  julia oszicar2magmom.jl OSZICAR
  julia oszicar2magmom.jl OSZICAR --type M_int
  julia oszicar2magmom.jl OSZICAR --output magmom.txt
        """
    )

    @add_arg_table! s begin
        "oszicar_file"
            help = "Path to OSZICAR file"
            required = true
        "--type", "-t"
            help = "Magnetic moment type to use (MW_int or M_int) [default: MW_int]"
            arg_type = String
            default = "MW_int"
            range_tester = (x -> x in ["MW_int", "M_int"])
        "--output", "-o"
            help = "Output filename (if not specified, output to stdout)"
            arg_type = String
            default = nothing
    end

    args = parse_args(s)

    # Check if file exists
    if !isfile(args["oszicar_file"])
        println(stderr, "Error: File not found: $(args["oszicar_file"])")
        return 1
    end

    # Parse OSZICAR file
    println(stderr, "Loading OSZICAR file: $(args["oszicar_file"])")
    println(stderr, "Magnetic moment type: $(args["type"])")

    data = parse_oszicar_last_step(args["oszicar_file"], args["type"])

    if isempty(data)
        println(stderr, "Error: Failed to extract data.")
        return 1
    end

    # Get the last step with data
    steps = sort(collect(keys(data)))
    if isempty(steps)
        println(stderr, "Error: No SCF steps found.")
        return 1
    end

    # Find the last step that has data
    last_step = nothing
    last_magmom = nothing
    for step in reverse(steps)
        if !isempty(data[step])
            last_step = step
            last_magmom = data[step]
            break
        end
    end

    if last_step === nothing || last_magmom === nothing
        println(stderr, "Error: No data found in any SCF step.")
        return 1
    end

    println(stderr, "Using last SCF step with data: $last_step")
    println(stderr, "Number of atoms: $(length(last_magmom) ÷ 3)")

    # Format and output
    magmom_line = format_magmom_line(last_magmom)

    if args["output"] !== nothing
        open(args["output"], "w") do f
            println(f, magmom_line)
        end
        println(stderr, "Output written to $(args["output"])")
    else
        println(magmom_line)
    end

    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

