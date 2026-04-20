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