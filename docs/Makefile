# Makefile for Magesty.jl documentation

.PHONY: build serve clean help

# Build documentation
build:
	julia --project=. make.jl

# Serve documentation locally (requires LiveServer.jl)
serve:
	julia --project=. -e "using LiveServer; serve(dir=\"build\")"

# Clean build directory
clean:
	rm -rf build/

# Open documentation in browser
open:
	open build/index.html

# Build and open
build-open: build open

# Help
help:
	@echo "Available targets:"
	@echo "  build       - Build documentation"
	@echo "  serve       - Serve documentation locally"
	@echo "  clean       - Clean build directory"
	@echo "  open        - Open documentation in browser"
	@echo "  build-open  - Build and open documentation"
	@echo "  help        - Show this help message"
