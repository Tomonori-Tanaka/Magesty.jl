# Installation

## Installing the Package

Open the Julia REPL and run:

```julia
using Pkg
Pkg.add(url="https://github.com/Tomonori-Tanaka/Magesty.jl")
```

## CLI Tools

Magesty.jl ships helper scripts under `tools/` (for example
`tools/vasp/vasp2extxyz.jl`). Run them directly with Julia against this
package's environment:

```sh
julia --project=@v$(VERSION.major).$(VERSION.minor) /path/to/Magesty.jl/tools/vasp/vasp2extxyz.jl ARGS
```

If you use these scripts frequently, wrap them in a shell function or
alias yourself.
