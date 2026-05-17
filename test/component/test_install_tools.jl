using Test
using Magesty

@testset "install_tools" begin
    mktempdir() do bindir
        # Suppress the "Installed: ..." / PATH hint stdout chatter; the
        # function returns nothing, so we only need to check side effects.
        redirect_stdout(devnull) do
            Magesty.install_tools(bindir = bindir)
        end

        expected = ["vasp2extxyz", "vasp2extxyz_recursive"]
        for name in expected
            wrapper = joinpath(bindir, name)
            @test isfile(wrapper)
            # Executable bit set for owner / group / other (0o755).
            @test (stat(wrapper).mode & 0o111) != 0
            # Wrapper is a thin shell script that delegates to julia.
            content = read(wrapper, String)
            @test startswith(content, "#!/bin/sh")
            @test occursin("exec julia", content)
        end
    end
end
