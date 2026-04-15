# Automatically installs CLI wrapper scripts into ~/.julia/bin/ on Pkg.add / Pkg.update.
pkg_dir = dirname(@__DIR__)

# Resolve dependencies into the package's own environment so that
# wrappers launched with --project="$pkg_dir" can find EzXML etc.
@info "Instantiating Magesty environment…"
run(`$(Base.julia_cmd()) --startup-file=no --project=$pkg_dir -e "import Pkg; Pkg.instantiate()"`)
@info "Done."

bindir  = joinpath(homedir(), ".julia", "bin")
mkpath(bindir)

tools = [
    "vasp2extxyz"           => joinpath("tools", "vasp", "vasp2extxyz.jl"),
    "vasp2extxyz_recursive" => joinpath("tools", "vasp", "vasp2extxyz_recursive.jl"),
]

for (name, rel_path) in tools
    script  = joinpath(pkg_dir, rel_path)
    wrapper = joinpath(bindir, name)
    write(wrapper, "#!/bin/sh\nexec julia --project=\"$pkg_dir\" \"$script\" \"\$@\"\n")
    chmod(wrapper, 0o755)
    @info "Installed: $wrapper"
end

@info "Add ~/.julia/bin to PATH if not already done:\n  export PATH=\"\$HOME/.julia/bin:\$PATH\""
