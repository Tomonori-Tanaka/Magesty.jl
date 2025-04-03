using Printf
function set_magmom_noncollinear(
	file::AbstractString,
	magmom::AbstractVector{<:Real},
	m_const::Bool = true,
)::String
	# Transform magmom to string form
	values_str = ""
	for (i, m_value) in enumerate(magmom)
        if i % 30 == 0
			values_str *= @sprintf("%.4f \\ \n", m_value)
		elseif i % 3 == 0
			values_str *= @sprintf("%.4f  ", m_value)
		else
			values_str *= @sprintf("%.4f ", m_value)
		end
	end

	magmom_str = "MAGMOM = $values_str"

    # Read file content
    content = read(file, String)
    
    # Split content into lines
    lines = split(content, '\n')
    new_lines = String[]
    
    i = 1
    while i <= length(lines)
        line = lines[i]
        # Skip lines containing MAGMOM
        if occursin("MAGMOM", line)
            # Check if the line contains backslash
            if i < length(lines) && occursin("\\", line)
                i += 1  # Skip next line
            end
        elseif occursin("M_CONSTR", line)
            if i < length(lines) && occursin("\\", line)
                i += 1  # Skip next line
            end
        else
            push!(new_lines, line)
        end
        i += 1
    end
    
    # Join lines to create new content
    new_content = join(new_lines, '\n')
    
    # Add MAGMOM and M_CONSTR if needed
    new_content = string(new_content, '\n', magmom_str)
    if m_const
        new_content = string(new_content, '\n', "M_CONSTR = $values_str")
    end
    
    return new_content
end


