test-all:
	julia --project -e 'using Pkg; Pkg.test()'

test-unit:
	TEST_MODE=unit julia --project -e 'using Pkg; Pkg.test()'

test-integration:
	TEST_MODE=integration julia --project -e 'using Pkg; Pkg.test()'

test-develop:
	TEST_MODE=develop julia --project -e 'using Pkg; Pkg.test()'

test-tools:
	julia tools/test/runtests.jl

test-jet:
	TEST_MODE=jet julia --project -e 'using Pkg; Pkg.test()'

test-aqua:
	TEST_MODE=aqua julia --project -e 'using Pkg; Pkg.test()'

bench-setup:
	julia --project=bench -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'

bench-sphericart:
	julia --project=bench bench/benchmark_sphericart.jl

bench-salcbasis:
	julia --project=bench bench/benchmark_salcbasis_hotspots.jl

bench-spherical-harmonics:
	julia --project=bench bench/benchmark_spherical_harmonics.jl

bench-threads:
	bash bench/run_threads_scaling.sh

test-sphericart:
	TEST_MODE=sphericart julia --project -e 'using Pkg; Pkg.test()'