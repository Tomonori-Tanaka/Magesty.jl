# Installation

## Installing the Package

Open the Julia REPL and run:

```julia
using Pkg
Pkg.add(url="https://github.com/Tomonori-Tanaka/Magesty.jl")
```

## Command-Line Interface

The `magesty` command-line tool is provided by the `MagestyCLI` package,
which lives in the `cli/` subdirectory of the repository. Because
`MagestyCLI` depends on the core `Magesty` package by path, it must be
built from a checkout of the repository, where the two packages sit side
by side. Clone the repository and build the launcher from the `cli/`
environment:

```sh
git clone https://github.com/Tomonori-Tanaka/Magesty.jl
cd Magesty.jl
julia --project=cli -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=cli cli/deps/build.jl
```

The first `julia` command makes the core `Magesty` package resolvable
inside the `cli/` environment; the second runs the build step that writes
the `magesty` launcher into `~/.julia/bin`. From a checkout, `make
install-cli` runs both steps in one command.

Add `~/.julia/bin` to your `PATH`:

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

The launcher runs in its own dedicated environment. After updating your
checkout (`git pull`), re-run the build step to refresh that environment:

```sh
julia --project=cli cli/deps/build.jl
```
