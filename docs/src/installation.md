# Installation

## Installing the Package

Open the Julia REPL and run:

```julia
using Pkg
Pkg.add(url="https://github.com/Tomonori-Tanaka/Magesty.jl")
```

## Command-Line Interface

The `magesty` command-line tool is provided by the `MagestyCLI` package,
which lives in the `cli/` subdirectory of the repository. Add it after
the core package, then build it:

```julia
using Pkg
Pkg.add(url="https://github.com/Tomonori-Tanaka/Magesty.jl")
Pkg.add(url="https://github.com/Tomonori-Tanaka/Magesty.jl", subdir="cli")
Pkg.build("MagestyCLI")
```

`Pkg.build("MagestyCLI")` writes the `magesty` launcher into
`~/.julia/bin`. Add that directory to your `PATH`:

```sh
export PATH="$HOME/.julia/bin:$PATH"
```

Then list the available commands:

```sh
magesty --help
```

See the [Tools](tools.md) page for the full list of `magesty` subcommands
and their options.

For example, convert a VASP run to extended XYZ:

```sh
magesty vasp extxyz vasprun.xml --oszicar OSZICAR --output frame.extxyz
```

The launcher runs in its own dedicated environment. After updating
Magesty.jl, re-run `Pkg.build("MagestyCLI")` to refresh that environment.
