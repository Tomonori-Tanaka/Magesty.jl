"""
    version.jl

Module for managing version information.
"""

module Version

export version_string

# Major version
const MAJOR = 0

# Minor version
const MINOR = 1

# Patch version
const PATCH = 0

# Prerelease information (set if needed)
const PRERELEASE = nothing

# Build metadata (set if needed)
const BUILD = nothing

"""
    version_string()

Returns the version information as a string.
"""
function version_string()
    v = "$MAJOR.$MINOR.$PATCH"
    if !isnothing(PRERELEASE)
        v *= "-$PRERELEASE"
    end
    if !isnothing(BUILD)
        v *= "+$BUILD"
    end
    return v
end

end # module 