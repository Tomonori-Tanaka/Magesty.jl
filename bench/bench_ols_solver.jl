using Magesty
using BenchmarkTools
using LinearAlgebra
using Random
using Printf

println("BLAS threads = ", BLAS.get_num_threads())
println("Julia threads = ", Threads.nthreads())

const SOLVE = Magesty.Fitting.solve_coefficients

function run_case(m::Int, p::Int; seed::Int = 1)
    rng = MersenneTwister(seed)
    X = randn(rng, m, p)
    y = randn(rng, m)

    b_ols   = SOLVE(OLS(), X, y)
    b_ridge = SOLVE(Ridge(lambda = 1e-8), X, y)

    res_ols   = norm(X * b_ols   - y)
    res_ridge = norm(X * b_ridge - y)

    t_ols = @belapsed $SOLVE(OLS(),                 $X, $y) samples=3 evals=1
    t_rdg = @belapsed $SOLVE(Ridge(lambda = 1e-8),  $X, $y) samples=3 evals=1

    @printf("m=%6d p=%4d  OLS=%7.3f s  Ridge=%7.3f s  ratio=%5.2fx  ‖res‖(OLS)=%.3e  ‖res‖(Ridge)=%.3e\n",
            m, p, t_ols, t_rdg, t_ols / t_rdg, res_ols, res_ridge)
end

for (m, p) in [(2000, 500), (5000, 1000), (10000, 1500)]
    run_case(m, p)
end
