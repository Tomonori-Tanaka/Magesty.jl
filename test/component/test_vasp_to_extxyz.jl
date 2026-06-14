"""
Regression tests for the `vasp_to_extxyz` API. The output must be
byte-identical to the golden extxyz fixtures, which were generated from
the VASP fixture inputs. The `magesty vasp extxyz` CLI subcommand is
covered separately by the `MagestyCLI` package tests.
"""

using Magesty
using Test

@testset "vasp_to_extxyz: $(system)" for system in ("FeRh", "IrMn3")
    dir          = joinpath(@__DIR__, "fixtures", system)
    vasprun_path = joinpath(dir, "vasprun.xml")
    oszicar_path = joinpath(dir, "OSZICAR")
    golden_text  = read(joinpath(dir, "$(system).extxyz"), String)

    # API function: the returned text is byte-identical to the golden file.
    @test vasp_to_extxyz(vasprun_path; oszicar = oszicar_path) == golden_text

    # The `output` keyword writes the same bytes to disk.
    api_outfile = tempname() * ".extxyz"
    try
        vasp_to_extxyz(vasprun_path; oszicar = oszicar_path, output = api_outfile)
        @test read(api_outfile, String) == golden_text
    finally
        isfile(api_outfile) && rm(api_outfile)
    end
end

# A non-converged run (scstep count >= NELM) is still written, tagged
# `converged=F`, and warns. We synthesize one by lowering NELM below the
# fixture's scstep count, so the count rule reports non-convergence.
@testset "vasp_to_extxyz: non-converged run warns and tags F" begin
    src = read(joinpath(@__DIR__, "fixtures", "FeRh", "vasprun.xml"), String)
    # FeRh runs 35 electronic steps; force the cap to 30 so 35 >= 30 → F.
    edited = replace(src, "name=\"NELM\">   100" => "name=\"NELM\">    30")
    @test edited != src   # guard: the substitution actually fired

    tmp = tempname() * ".xml"
    try
        write(tmp, edited)
        text = @test_logs (:warn, r"did not converge") vasp_to_extxyz(tmp)
        header = split(text, "\n")[2]
        @test occursin("converged=F", header)
    finally
        isfile(tmp) && rm(tmp)
    end
end

# A parse failure must name the offending file path, so a batch conversion
# points the user straight at the bad input.
@testset "vasp_to_extxyz: parse error names the file path" begin
    # Missing file: the path appears in the error message.
    missing_path = joinpath(tempdir(), "does_not_exist_$(getpid()).xml")
    @test !isfile(missing_path)   # guard
    err = try
        vasp_to_extxyz(missing_path)
        nothing
    catch e
        e
    end
    @test err isa Magesty.VaspParseError
    @test occursin(missing_path, sprint(showerror, err))

    # Malformed XML mid-parse: still a VaspParseError carrying the path.
    tmp = tempname() * ".xml"
    try
        write(tmp, "<modeling><not-valid-vasprun></modeling>")
        err2 = try
            vasp_to_extxyz(tmp)
            nothing
        catch e
            e
        end
        @test err2 isa Magesty.VaspParseError
        @test occursin(tmp, sprint(showerror, err2))
    finally
        isfile(tmp) && rm(tmp)
    end

    # The OSZICAR path is named when magnetic parsing fails. The vasprun is
    # valid; the OSZICAR is missing, so the OSZICAR path must surface.
    valid_vasprun = joinpath(@__DIR__, "fixtures", "FeRh", "vasprun.xml")
    bad_oszicar   = joinpath(tempdir(), "no_such_OSZICAR_$(getpid())")
    @test !isfile(bad_oszicar)   # guard
    err3 = try
        vasp_to_extxyz(valid_vasprun; oszicar = bad_oszicar)
        nothing
    catch e
        e
    end
    @test err3 isa Magesty.VaspParseError
    @test occursin(bad_oszicar, sprint(showerror, err3))
end

println("All vasp_to_extxyz tests passed.")
