"""
	Optimize.jl

This module contains functions for optimizing the SCE coefficients.
"""

module Optimize
using ..BasisSet

export optimize_sce

struct SCEOptimizer
    SCE::Vector{Float64}
end

end
