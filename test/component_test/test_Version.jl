using Test
using TOML
using Magesty

# Guard against the version constants drifting away from Project.toml
# (the historical reason for adding this module to the refactor sweep:
# `src/common/version.jl` used to hard-code MAJOR/MINOR/PATCH/PRERELEASE
# next to Project.toml's own `version` field, and the two could fall out
# of sync silently).

@testset "Version" begin
    project_toml_path = joinpath(pkgdir(Magesty), "Project.toml")
    project_version = TOML.parsefile(project_toml_path)["version"]::String

    @test Magesty.Version.version_string() == project_version
    # SemVer shape: MAJOR.MINOR.PATCH with optional -PRERELEASE / +BUILD tags.
    @test occursin(
        r"^[0-9]+\.[0-9]+\.[0-9]+(-[A-Za-z0-9.]+)?(\+[A-Za-z0-9.]+)?$",
        Magesty.Version.version_string(),
    )
end
