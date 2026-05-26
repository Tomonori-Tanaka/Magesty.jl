using Test
using Magesty
import TOML

# Fixtures -------------------------------------------------------------

const DIMER_TOML_IO = joinpath(@__DIR__, "..", "integration", "dimer", "input.toml")

# FM / AFM spin configurations for the 2-atom dimer. The component values
# (0.0, ±1.0, ±0.1) are exact in binary floating point, so the EMBSET
# round trip in "signature families agree" recovers them bit-for-bit; do
# not replace them with arbitrary floats or that byte-identity check breaks.
function _dimer_configs_io()
    SC = Magesty.SpinConfigs.SpinConfig
    fm = SC(-1.0, [1.0, 1.0],
            [0.0 0.0; 0.0 0.0; 1.0 1.0],
            [0.1 0.0; 0.0 -0.1; 0.0 0.0])
    afm = SC(1.0, [1.0, 1.0],
             [0.0 0.0; 0.0 0.0; -1.0 1.0],
             [-0.1 0.0; 0.0 0.1; 0.0 0.0])
    return [fm, afm]
end

# Write a minimal EMBSET file for the dimer configs. Each config block is a
# comment line, an energy line, and one line per atom holding the spin
# moment vector (magnitude x direction) followed by the local field.
function _write_dimer_embset(path::AbstractString, configs)
    open(path, "w") do io
        for (c, sc) in enumerate(configs)
            println(io, "# $c, test, energy unit = eV")
            println(io, sc.energy)
            for a = 1:length(sc.magmom_size)
                m = sc.magmom_size[a] .* sc.spin_directions[:, a]
                b = sc.local_magfield[:, a]
                println(io, "  $a  $(m[1]) $(m[2]) $(m[3])  $(b[1]) $(b[2]) $(b[3])")
            end
        end
    end
    return path
end

# Non-comment, non-empty lines of a written file.
_data_lines(text::AbstractString) =
    [ln for ln in split(text, '\n')
     if !isempty(strip(ln)) && !startswith(strip(ln), "#")]

# Tests ----------------------------------------------------------------

@testset "FitCheckIO" begin
    input = TOML.parsefile(DIMER_TOML_IO)
    basis = SCEBasis(input; verbosity = false)
    configs = _dimer_configs_io()
    dataset = SCEDataset(basis, configs)
    # `Ridge(lambda = 1e-8)` rather than `OLS()` because the dimer
    # FM / AFM configs have spins along z, which makes every torque
    # zero by S × S = 0; at the default `torque_weight = 1.0` the
    # design matrix is then entirely zero. The test only checks file
    # format and the round-trip of predict_* outputs, so the
    # regularizer is incidental.
    f = fit(SCEFit, dataset, Ridge(lambda = 1e-8); verbosity = false)
    n_atoms = basis.structure.supercell.num_atoms

    @testset "energy file format" begin
        path = joinpath(mktempdir(), "energy_list.txt")
        @test write_energies(f, path) === nothing
        text = read(path, String)
        @test occursin("# data index", text)
        @test occursin("# unit of energy is eV", text)

        rows = _data_lines(text)
        @test length(rows) == length(configs)
        pred = predict_energy(f, configs)
        for (i, row) in enumerate(rows)
            cols = split(row)
            @test length(cols) == 3                       # index, DFT, SCE
            @test parse(Int, cols[1]) == i
            @test parse(Float64, cols[2]) == configs[i].energy
            @test parse(Float64, cols[3]) ≈ pred[i] atol = 1e-9
        end
    end

    @testset "torque file format" begin
        path = joinpath(mktempdir(), "torque_list.txt")
        @test write_torques(f, path) === nothing
        text = read(path, String)
        @test occursin("# atom index", text)
        @test occursin("# unit of torque is eV", text)
        @test count(ln -> startswith(strip(ln), "# data index:"),
                    split(text, '\n')) == length(configs)

        rows = _data_lines(text)
        # One row per atom per config (the x/y/z components are columns).
        @test length(rows) == length(configs) * n_atoms
        @test all(length(split(row)) == 8 for row in rows)

        # First block, first atom: columns match observed / predicted.
        pred = predict_torque(f, configs)
        cols = split(rows[1])
        @test parse(Int, cols[1]) == 1
        for k = 1:3
            @test parse(Float64, cols[2+k]) ≈ configs[1].torques[k, 1] atol = 1e-12
            @test parse(Float64, cols[5+k]) ≈ pred[1][k, 1] atol = 1e-9
        end
    end

    @testset "default filenames" begin
        cd(mktempdir()) do
            write_energies(f)
            write_torques(f)
            @test isfile("energy_list.txt")
            @test isfile("torque_list.txt")
        end
    end

    @testset "signature families agree" begin
        dir = mktempdir()
        embset = _write_dimer_embset(joinpath(dir, "EMBSET"), configs)
        model = SCEModel(f)

        # Energy: the SCEFit-self form and the (predictor, data) form with
        # each accepted data type must produce byte-identical files.
        e_fit = joinpath(dir, "e_fit.txt"); write_energies(f, e_fit)
        e_ds  = joinpath(dir, "e_ds.txt");  write_energies(model, dataset, e_ds)
        e_cfg = joinpath(dir, "e_cfg.txt"); write_energies(model, configs, e_cfg)
        e_emb = joinpath(dir, "e_emb.txt"); write_energies(model, embset, e_emb)
        e_fitp = joinpath(dir, "e_fitp.txt"); write_energies(f, configs, e_fitp)
        ref_e = read(e_fit, String)
        @test read(e_ds, String) == ref_e
        @test read(e_cfg, String) == ref_e
        @test read(e_emb, String) == ref_e
        @test read(e_fitp, String) == ref_e          # SCEFit as predictor

        # Torque: same check.
        t_fit = joinpath(dir, "t_fit.txt"); write_torques(f, t_fit)
        t_ds  = joinpath(dir, "t_ds.txt");  write_torques(model, dataset, t_ds)
        t_cfg = joinpath(dir, "t_cfg.txt"); write_torques(model, configs, t_cfg)
        t_emb = joinpath(dir, "t_emb.txt"); write_torques(model, embset, t_emb)
        t_fitp = joinpath(dir, "t_fitp.txt"); write_torques(f, configs, t_fitp)
        ref_t = read(t_fit, String)
        @test read(t_ds, String) == ref_t
        @test read(t_cfg, String) == ref_t
        @test read(t_emb, String) == ref_t
        @test read(t_fitp, String) == ref_t          # SCEFit as predictor
    end
end
