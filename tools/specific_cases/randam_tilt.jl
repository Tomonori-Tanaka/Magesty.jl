include("../vasptools.jl")

using .VaspTools
using Random
using ArgParse
using LinearAlgebra

"""
    random_tilt(
        incar::String,
        num_atoms::Integer,
        theta_range::Tuple{Float64, Float64},
        num_patterns::Integer,
        ;
        file_prefix::String = "pattern",
    )

Generate random tilt MAGMOM and M_CONSTR patterns for a VASP calculation.

# Arguments
- `incar::String`: Path to the input INCAR file
- `num_atoms::Integer`: Number of atoms in the structure
- `theta_range::Tuple{Float64, Float64}`: Range of tilt angles in degrees from the original magnetic moment direction (min, max)
- `num_patterns::Integer`: Number of random patterns to generate
- `file_prefix::String`: Prefix for output files (default: "pattern")

# Returns
- Generates `num_patterns` INCAR files with random magnetic moment orientations
- Each file is named as `{file_prefix}_{pattern_number}.incar`
- Pattern numbers are zero-padded based on the total number of patterns

# Notes
- The tilt angle θ is defined as the angle from the original magnetic moment direction
- The azimuthal angle φ is defined as the angle around the original magnetic moment direction
- The magnitude of each magnetic moment is preserved
- The rotation axis is randomly chosen perpendicular to the original magnetic moment
- Atoms with zero magnetic moment are kept unchanged
- The original magnetic moment direction is taken from M_CONSTR in the input INCAR file
- Both MAGMOM and M_CONSTR in the output files are set to the same values
"""
function random_tilt(
    incar::String,
    num_atoms::Integer,
    theta_range::Tuple{Float64, Float64},
    num_patterns::Integer,
    ;
    file_prefix::String = "pattern",
)
    # Parse INCAR file
    incar_dict = parse_incar(incar)
    
    # Get M_CONSTR values
    m_constr = incar_dict[:M_CONSTR]
    
    # Calculate required number of digits for zero-padding
    digits = length(string(num_patterns))
    
    # Convert theta_range from degrees to radians
    theta_range_rad = (theta_range[1] * π/180, theta_range[2] * π/180)
    
    # Generate random patterns
    for i in 1:num_patterns
        new_magmom = Vector{Float64}()
        new_m_constr = Vector{Float64}()
        sizehint!(new_magmom, 3 * num_atoms)
        sizehint!(new_m_constr, 3 * num_atoms)
        
        # Generate random orientations for each atom
        for j in 1:num_atoms
            # Get original magnetic moment components
            m_orig = [m_constr[3*(j-1)+1], m_constr[3*(j-1)+2], m_constr[3*(j-1)+3]]
            
            # Calculate magnitude of original magnetic moment
            m = norm(m_orig)
            
            if m ≈ 0.0
                # Keep zero magnetic moment atoms unchanged
                append!(new_magmom, m_orig)
                append!(new_m_constr, m_orig)
            else
                # Normalize original magnetic moment
                m_hat = m_orig / m
                
                # Generate random tilt angle and azimuthal angle
                theta = rand() * (theta_range_rad[2] - theta_range_rad[1]) + theta_range_rad[1]
                phi = rand() * 2π
                
                # Find a vector perpendicular to m_hat
                if abs(m_hat[1]) < 0.5
                    v = [1.0, 0.0, 0.0]
                else
                    v = [0.0, 1.0, 0.0]
                end
                
                # Make it perpendicular to m_hat using Gram-Schmidt
                axis1 = v - dot(v, m_hat) * m_hat
                axis1 = axis1 / norm(axis1)
                
                # Create a second perpendicular axis
                axis2 = cross(m_hat, axis1)
                
                # Create rotation axis using azimuthal angle
                axis = cos(phi) * axis1 + sin(phi) * axis2
                
                # Rodrigues' rotation formula
                cos_theta = cos(theta)
                sin_theta = sin(theta)
                m_new = cos_theta * m_orig + sin_theta * cross(axis, m_orig) + (1 - cos_theta) * dot(axis, m_orig) * axis
                
                # Append new values
                append!(new_magmom, m_new)
                append!(new_m_constr, m_new)
            end
        end
        
        # Create new INCAR
        new_incar_dict = copy(incar_dict)
        new_incar_dict[:MAGMOM] = new_magmom
        new_incar_dict[:M_CONSTR] = new_m_constr
        
        # Write to file with zero-padded pattern number
        output_file = "$(file_prefix)_$(lpad(i, digits, "0")).incar"
        write_incar(output_file, new_incar_dict)
    end
end

function parse_commandline()
    s = ArgParseSettings(
        description = "Generate random tilt MAGMOM and M_CONSTR patterns for VASP calculations.",
        version = "1.0.0",
        add_version = true,
    )

    @add_arg_table! s begin
        "--incar", "-i"
            help = "Path to the input INCAR file"
            required = true
            arg_type = String

        "--num-atoms", "-n"
            help = "Number of atoms in the structure"
            required = true
            arg_type = Int

        "--theta-min", "-t"
            help = "Minimum tilt angle in degrees from the original magnetic moment direction"
            required = true
            arg_type = Float64

        "--theta-max", "-T"
            help = "Maximum tilt angle in degrees from the original magnetic moment direction"
            required = true
            arg_type = Float64

        "--num-patterns", "-p"
            help = "Number of random patterns to generate"
            required = true
            arg_type = Int

        "--prefix", "-f"
            help = "Prefix for output files"
            default = "pattern"
            arg_type = String
    end

    return parse_args(s)
end

function main()
    # Parse command line arguments
    args = parse_commandline()
    
    # Call random_tilt function with parsed arguments
    random_tilt(
        args["incar"],
        args["num-atoms"],
        (args["theta-min"], args["theta-max"]),
        args["num-patterns"],
        file_prefix = args["prefix"],
    )
end

# Execute main function if script is run directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end