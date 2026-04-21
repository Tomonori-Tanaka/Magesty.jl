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

bench-sphericart:
	TEST_MODE=bench_sphericart julia --project -e 'using Pkg; Pkg.test()'

test-sphericart:
	TEST_MODE=sphericart julia --project -e 'using Pkg; Pkg.test()'

bench-optimize-sphericart:
	TEST_MODE=bench_optimize_sphericart julia --project -e 'using Pkg; Pkg.test()'

bench-optimize-sphericart-fept:
	TEST_MODE=bench_optimize_sphericart BENCH_EXAMPLE=fept_tetragonal_2x2x2 julia --project -e 'using Pkg; Pkg.test()'