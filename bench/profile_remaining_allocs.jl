#!/usr/bin/env julia
# Per-call allocation count for the Fitting design-matrix inner loops.
# Wraps `design_matrix_energy_element` and the cluster-major torque
# accumulator and counts heap allocations under repeated invocation,
# separating the per-call overhead from one-time setup. Useful when
# investigating a regression in scratch-buffer reuse.

using Printf
using Magesty
using Magesty.Fitting: design_matrix_energy_element,
                       _accumulate_grad_torque_cluster!,
                       build_sh_cache_energy, build_sh_cache_torque,
                       _compute_l_max
using Magesty.SpinConfigs: read_embset
using StaticArrays

const EXAMPLE_DIR = joinpath(@__DIR__, "..", "test", "integration", "fept_tetragonal_2x2x2")

basis = Magesty.load(SCEBasis, joinpath(EXAMPLE_DIR, "system.xml"))
spinconfigs = read_embset(joinpath(EXAMPLE_DIR, "EMBSET"))

l_max = _compute_l_max(basis.salcbasis.salc_list)
sh_cache_e = build_sh_cache_energy(spinconfigs[1].spin_directions, l_max)
sh_cache_t = build_sh_cache_torque(spinconfigs[1].spin_directions, l_max)
num_atoms = size(spinconfigs[1].spin_directions, 2)
num_salcs = length(basis.salcbasis.salc_list)

# Warm up
for cbc in basis.salcbasis.salc_list[1]
    design_matrix_energy_element(cbc, sh_cache_e)
end

# Per-call alloc count for energy_element
function bench_one_call(cbc, sh_cache, n)
    allocs = @allocations begin
        for _ in 1:n
            design_matrix_energy_element(cbc, sh_cache)
        end
    end
    return allocs / n
end

println("=== Per-call allocations for design_matrix_energy_element (averaged over n=1000) ===")
for (salc_idx, group) in enumerate(basis.salcbasis.salc_list[1:min(5, end)])
    for (cbc_idx, cbc) in enumerate(group)
        N = length(cbc.atoms)
        a = bench_one_call(cbc, sh_cache_e, 1000)
        @printf "  salc=%-2d cbc=%-2d N=%d (R=%d) ls=%-8s → %.1f allocs/call\n" salc_idx cbc_idx N (N+1) string(cbc.ls) a
    end
end

# Per-call for _accumulate_grad_torque_cluster!
println("\n=== Per-call allocations for _accumulate_grad_torque_cluster! ===")
grad_buf = zeros(Float64, 3, num_atoms, num_salcs)
for (salc_idx, group) in enumerate(basis.salcbasis.salc_list[1:min(5, end)])
    for (cbc_idx, cbc) in enumerate(group)
        # `@inbounds` here matches the production call context inside
        # `build_design_matrix_torque`, so the `@boundscheck @assert
        # allunique(cluster_atoms)` inside the kernel is elided and the
        # reported alloc count reflects production behavior.
        @inbounds for _ in 1:50  # warm
            _accumulate_grad_torque_cluster!(grad_buf, cbc, salc_idx, sh_cache_t)
        end
        allocs = @allocations begin
            @inbounds for _ in 1:1000
                _accumulate_grad_torque_cluster!(grad_buf, cbc, salc_idx, sh_cache_t)
            end
        end
        a = allocs / 1000
        N = length(cbc.atoms)
        @printf "  salc=%-2d cbc=%-2d N=%d ls=%-8s → %.1f allocs/call\n" salc_idx cbc_idx N string(cbc.ls) a
    end
end

using Printf
