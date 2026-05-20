# Command-line interface for Magesty, built with Comonicon.
#
# This file is included last in `module Magesty` so the `@cast` / `@main`
# macros see the full package API. They define `Magesty.command_main`
# (the `magesty` entry point) and `Magesty.comonicon_install` (writes the
# launcher into ~/.julia/bin); the latter is invoked from `deps/build.jl`.
#
# Subcommands are added as `@cast` functions; `vasp2extxyz` and the rest
# of the VASP converters land in later milestones.

# `@cast` / `@main` are qualified with `Comonicon.` because `@main`
# would otherwise clash with `Base.@main`.
import Comonicon

"""
Print the installed Magesty.jl version.
"""
Comonicon.@cast function version()
	println(pkgversion(@__MODULE__))
	return nothing
end

"""
Magesty command-line interface.

Build general effective spin models from noncollinear spin DFT
calculations using the spin-cluster expansion.
"""
Comonicon.@main
