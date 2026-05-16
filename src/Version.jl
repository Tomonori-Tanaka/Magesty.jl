"""
    version.jl

Module for accessing Magesty's version. The version itself is owned by
`Project.toml`; this module is a thin runtime accessor so that release
bumps require editing one file only.
"""
module Version

export version_string

"""
    version_string() -> String

Return the package version as a string by reading it from
`Project.toml` via `pkgversion`. Bumping the version is therefore a
single-file edit in `Project.toml`; no constants need to be kept in
sync here.

Returns `"unknown"` if the module was loaded outside of a package
context (e.g., included directly from a script).
"""
function version_string()::String
    v = pkgversion(@__MODULE__)
    return v === nothing ? "unknown" : string(v)
end

end # module
