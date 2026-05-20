# Installation

## Installing the Package

Open the Julia REPL and run:

```julia
using Pkg
Pkg.add(url="https://github.com/Tomonori-Tanaka/Magesty.jl")
```

## Command-Line Interface

Magesty.jl provides a `magesty` command-line tool. Install it after
adding the package:

```julia
using Pkg
Pkg.build("Magesty")
```

This writes the `magesty` launcher into `~/.julia/bin`. Add that
directory to your `PATH`:

```sh
export PATH="$HOME/.julia/bin:$PATH"
```

Then list the available commands:

```sh
magesty --help
```

For example, convert a VASP run to extended XYZ:

```sh
magesty vasp extxyz vasprun.xml --oszicar OSZICAR --output frame.extxyz
```

The launcher resolves Magesty.jl by name through its own environment, so
it keeps working after the package is updated — no reinstall needed.
