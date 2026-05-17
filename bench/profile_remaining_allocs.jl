#!/usr/bin/env julia
# Per-call allocation count for the Fitting design-matrix inner loops.
# Wraps `design_matrix_energy_element` and `calc_∇ₑu!` and counts
# heap allocations under repeated invocation, separating the per-call
# overhead from one-time setup. Useful when investigating a regression
# in scratch-buffer reuse.

using Printf
using Magesty
using Magesty.Fitting: design_matrix_energy_element, calc_∇ₑu!,
                       EnergyWorkspace, GradWorkspace
using Magesty.SpinConfigs: read_embset
using StaticArrays

const EXAMPLE_DIR = joinpath(@__DIR__, "..", "test", "integration", "fept_tetragonal_2x2x2")

basis = Magesty.load(SCEBasis, joinpath(EXAMPLE_DIR, "system.xml"))
spinconfigs = read_embset(joinpath(EXAMPLE_DIR, "EMBSET.dat"))

ws_e = EnergyWorkspace()
ws_g = GradWorkspace()

# Warm up
for cbc in basis.salcbasis.salc_list[1]
    design_matrix_energy_element(cbc, spinconfigs[1].spin_directions, basis.symmetry, ws_e)
end

# Per-call alloc count for energy_element
function bench_one_call(cbc, sd, sym, ws, n)
    allocs = @allocations begin
        for _ in 1:n
            design_matrix_energy_element(cbc, sd, sym, ws)
        end
    end
    return allocs / n
end

println("=== Per-call allocations for design_matrix_energy_element (averaged over n=1000) ===")
for (salc_idx, group) in enumerate(basis.salcbasis.salc_list[1:min(5, end)])
    for (cbc_idx, cbc) in enumerate(group)
        N = length(cbc.atoms)
        a = bench_one_call(cbc, spinconfigs[1].spin_directions, basis.symmetry, ws_e, 1000)
        @printf "  salc=%-2d cbc=%-2d N=%d (R=%d) ls=%-8s → %.1f allocs/call\n" salc_idx cbc_idx N (N+1) string(cbc.ls) a
    end
end

# Per-call for calc_∇ₑu!
println("\n=== Per-call allocations for calc_∇ₑu! ===")
result = MVector{3, Float64}(0.0, 0.0, 0.0)
for (salc_idx, group) in enumerate(basis.salcbasis.salc_list[1:min(5, end)])
    for (cbc_idx, cbc) in enumerate(group)
        for _ in 1:50  # warm
            calc_∇ₑu!(result, cbc, 1, spinconfigs[1].spin_directions, basis.symmetry, ws_g)
        end
        allocs = @allocations begin
            for _ in 1:1000
                calc_∇ₑu!(result, cbc, 1, spinconfigs[1].spin_directions, basis.symmetry, ws_g)
            end
        end
        a = allocs / 1000
        N = length(cbc.atoms)
        @printf "  salc=%-2d cbc=%-2d N=%d ls=%-8s → %.1f allocs/call\n" salc_idx cbc_idx N string(cbc.ls) a
    end
end

using Printf
