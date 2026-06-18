# Installation

## Installing the Package

Magesty is registered in the Julia General registry. Open the Julia REPL
and run:

```julia
using Pkg
Pkg.add("Magesty")
```

To track the latest development version instead, add it by URL:

```julia
using Pkg
Pkg.add(url="https://github.com/Tomonori-Tanaka/Magesty.jl")
```

## Command-Line Interface

The `magesty` command-line tool is provided by the `MagestyCLI` package in
the `cli/` subdirectory of the repository. Installing it writes a small
launcher script into `~/.julia/bin`; the launcher runs in its own
dedicated environment, independent of your active project.

### Install without cloning

Because the core `Magesty` package is in the General registry, the CLI can
be installed directly from the repository's `cli/` subdirectory. Run the
command in a temporary environment so it does not touch your default
project:

```julia
using Pkg
Pkg.activate(temp = true)
Pkg.add(url = "https://github.com/Tomonori-Tanaka/Magesty.jl", subdir = "cli")
```

Adding the package automatically runs the build step that writes the
`magesty` launcher into `~/.julia/bin`. The launcher is self-contained, so
the temporary environment can be discarded afterwards. To upgrade later,
re-run the same two commands; they fetch the latest revision and rebuild
the launcher.

### Install from a checkout (development)

When developing Magesty, build the CLI against a local checkout so the
launcher tracks your working copy:

```sh
git clone https://github.com/Tomonori-Tanaka/Magesty.jl
cd Magesty.jl
julia --project=cli -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=cli cli/deps/build.jl
```

The first `julia` command makes the core `Magesty` package resolvable
inside the `cli/` environment by path; the second runs the build step that
writes the launcher. `make install-cli` runs both steps in one command.
After updating your checkout (`git pull`), re-run the build step to refresh
the launcher's environment:

```sh
julia --project=cli cli/deps/build.jl
```

### Using the launcher

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
