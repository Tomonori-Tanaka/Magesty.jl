using JET
using Test
using Magesty

@testset "JET static analysis" begin
    # `report_package` macro-expands every top-level form in the package
    # through JET's JuliaInterpreter-based loader. Comonicon's `@main`
    # macro (in `src/CLI.jl`) runs command-tree code generation at
    # expansion time that fails under that loader
    # (`Comonicon.AST.Entry has no field name`). This is an external
    # Comonicon/JET incompatibility, not a Magesty defect — `using
    # Magesty` and `make test-all` are unaffected. The fix is structural:
    # extracting the CLI into a separate package removes Comonicon from
    # the `Magesty` source tree that `report_package` walks. Until then
    # this test is marked broken; it returns to a real assertion
    # automatically once `report_package` stops throwing.
    result = try
        report_package(
            Magesty;
            toplevel_logger = nothing,
            target_modules = (Magesty,),
        )
    catch err
        @warn "JET could not analyze Magesty; skipping (Comonicon `@main` " *
              "is incompatible with JET's loader)" exception = err
        nothing
    end

    if result === nothing
        @test_broken false
    else
        reports = JET.get_reports(result)
        n = length(reports)

        if n > 0
            show(result)
            @warn "JET found $n potential issue(s) in Magesty"
        else
            @info "JET: no issues found"
        end

        @test n == 0
    end
end
