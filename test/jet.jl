using JET
using Test
using Magesty

@testset "JET static analysis" begin
    # Analyze Magesty's own code only
    result = report_package(
        Magesty;
        toplevel_logger = nothing,
        target_modules = (Magesty,),
    )

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
