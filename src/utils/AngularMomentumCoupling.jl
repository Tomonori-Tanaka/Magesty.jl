module AngularMomentumCoupling

# AngularMomentumCoupling.jl
# Build coupled angular momentum bases for arbitrary N bodies.
# - Enumerate all left-coupling paths and final total L (Lf)
# - Build coefficient tensors C^{(Lf,M)} by chaining Wigner 3j (complex basis)
# - Convert to real (tesseral) basis for sites and the final multiplet
#
# Dependencies:
#   ] add WignerSymbols
#
# Public API (exported):
#   mrange
#   enumerate_paths_left_all
#   coeff_tensor_complex
#   complex_to_real_tensor
#   build_all_complex_bases
#   build_all_real_bases
#
# Notes:
# - This implementation uses integer l only (no half-integers).
# - Left-coupling tree is used; other trees are related by Racah (6j) recoupling.
# - Inside a fixed tree, different (L12, L123, ..., Lf) labels produce
#   orthogonal coupled states (up to floating round-off).

using LinearAlgebra
using WignerSymbols  # provides clebschgordan
using OffsetArrays   # enables direct indexing by magnetic quantum numbers m

using ..SphericalHarmonicsTransforms

# ---------------- Utilities ----------------

"""
	mrange(l::Int) -> Vector{Int}

Return the integer m range for a given l: -l, ..., 0, ..., +l.
"""
mrange(l::Int) = collect((-l):l)


"""
	nmode_mul(A, M, n) -> Array

Mode-n product: apply matrix M on mode n (1-based) of N-way tensor A.
Equivalent to unfolding_n(C) = M * unfolding_n(A).
"""
function nmode_mul(A::AbstractArray, M::AbstractMatrix, n::Int)
	sz = size(A)
	d  = length(sz)
	@assert 1 <= n <= d "invalid mode index n=$n for tensor of order $d"

	# Permute so that mode n is first
	perm = (n, setdiff(1:d, (n,))...)
	invp = invperm(perm)

	B    = permutedims(A, perm)
	Bmat = reshape(B, sz[n], :)
	Cmat = M * Bmat
	C    = reshape(Cmat, (size(M, 1), sz[setdiff(1:d, (n,))]...))
	return permutedims(C, invp)
end

# ---------------- 1) Enumerate all left-coupling paths and final Lf ----------------

"""
	enumerate_paths_left_all(ls::Vector{Int}) -> Vector{Tuple{Vector{Int}, Int}}

Enumerate all admissible left-coupling sequences (triangle rules) for
ls = [l1, l2, ..., lN].

Return a vector of (Lseq, Lf):
  - Lseq :: Vector{Int} is [L12, L123, ..., L_{1...(N-1)}] WITHOUT the final stage
		   (length N-2 when N>=2; empty when N=2).
  - Lf   :: Int is the final total angular momentum after coupling (J_f).

Special cases:
  - N == 1: return [(Int[], ls[1])].
  - N == 2: return [(Int[], Lf)] for all Lf in |l1-l2| : (l1+l2).
"""
function enumerate_paths_left_all(ls::Vector{Int})::Vector{Tuple{Vector{Int}, Int}}
	N = length(ls)
	@assert N >= 1 "ls must be non-empty; got N=$N"

	paths = Tuple{Vector{Int}, Int}[]

	# N = 1: no coupling, final Lf must equal l1
	if N == 1
		push!(paths, (Int[], ls[1]))
		return paths
	end

	# N = 2: Lseq is empty, list all Lf allowed by triangle rule
	if N == 2
		for Lf in abs(ls[1]-ls[2]):(ls[1]+ls[2])
			push!(paths, (Int[], Lf))
		end
		return paths
	end

	# N >= 3: recurse through intermediate stages
	function recbuild(k::Int, Lprev::Int, seq::Vector{Int})
		if k == N-1
			# Final stage: (Lprev, lN) -> Lf for all triangle-allowed Lf
			for Lf in abs(Lprev-ls[end]):(Lprev+ls[end])
				push!(paths, (copy(seq), Lf))
			end
			return
		end
		lk1 = ls[k+1]
		for L in abs(Lprev-lk1):(Lprev+lk1)  # triangle rule
			recbuild(k+1, L, vcat(seq, L))
		end
	end

	# First stage: (l1, l2) -> L12
	for L12 in abs(ls[1]-ls[2]):(ls[1]+ls[2])
		recbuild(2, L12, [L12])
	end
	return paths
end

# ---------------- 2) Build complex coefficient tensor C^{(Lf,M)} ----------------

"""
	coeff_tensor_complex(ls, Lseq, Lf; T=Float64) -> Array{T, N+1}

Given ls = [l1,...,lN], an intermediate sequence Lseq = [L12, L123, ...] (length N-2),
and final Lf, build the (N+1)-way tensor C^(Lf,M) over m_1..m_N and final Mf (complex basis):

  C[i1, i2, ..., iN, q] = < (l1,m1) ... (lN,mN) | (Lf, Mf) >

where i_k enumerates m_k = -l_k..+l_k, and q enumerates Mf = -Lf..+Lf.

Construction uses chaining of Clebsch–Gordan coefficients.
"""
function coeff_tensor_complex(
	ls::Vector{Int},
	Lseq::Vector{Int},
	Lf::Int;
	T::Type{<:Real} = Float64,
)
	N = length(ls)
	mr = [mrange(l) for l in ls]
	dims = [2*l + 1 for l in ls]
	dims_all = vcat(dims, 2*Lf + 1)
	C = zeros(T, dims_all...)  # order N+1: m_1..m_N, Mf

	# N = 1: delta coupling
	if N == 1
		l1 = ls[1]
		m1 = mr[1]
		for (i, mi) in enumerate(m1), (q, Mf) in enumerate(mrange(Lf))
			C[i, q] = (Lf == l1 && Mf == mi) ? one(T) : zero(T)
		end
		return C
	end

	# N = 2: single CG coefficient
	if N == 2
		m1, m2 = mr[1], mr[2]
		for (i, mi) in enumerate(m1), (j, mj) in enumerate(m2), (q, Mf) in enumerate(mrange(Lf))
			C[i, j, q] = float(clebschgordan(ls[1], mi, ls[2], mj, Lf, Mf))
		end
		return C
	end

	# N >= 3: forward dynamic programming over intermediate M's
	Lprev = Lseq[end]
	inds  = CartesianIndices(Tuple(dims))  # Convert Vector to Tuple for multi-dimensional indexing
	for I in inds
		# Extract (m1,...,mN) for this multi-index
		ms = ntuple(i -> mr[i][I[i]], N)

		# Stage 1: (l1,l2) -> L12 (sum over M12)
		L12 = Lseq[1]
		acc_prev = Dict{Int, T}()    # map: M12 -> amplitude
		for M12 in (-L12):L12
			v = float(clebschgordan(ls[1], ms[1], ls[2], ms[2], L12, M12))
			if v != 0
				acc_prev[M12] = v
			end
		end

		# Intermediate stages: (L_{...}, l_{t+1}) -> L_t
		for t in 2:(N-2)
			Lt_1, Lt = Lseq[t-1], Lseq[t]
			mt1 = ms[t+1]
			acc = Dict{Int, T}()
			for (Mprev, aval) in acc_prev
				for Mt in (-Lt):Lt
					v = aval * float(clebschgordan(Lt_1, Mprev, ls[t+1], mt1, Lt, Mt))
					if v != 0
						acc[Mt] = get(acc, Mt, zero(T)) + v
					end
				end
			end
			acc_prev = acc
		end

		# Final stage: (Lprev, lN) -> Lf ; sum over Mprev
		mN = ms[end]
		for (q, Mf) in enumerate(mrange(Lf))
			total_amplitude = zero(T)
			for (Mprev, aval) in acc_prev
				total_amplitude += aval * float(clebschgordan(Lprev, Mprev, ls[end], mN, Lf, Mf))
			end
			C[I, q] = total_amplitude
		end
	end
	return C
end

# ---------------- 3) Complex -> Real (tesseral) conversion ----------------

"""
	complex_to_real_tensor(Ccx, ls, Lf) -> Array{Float64, N+1}

Apply unitary transforms (site-wise and final-Lf) to convert from complex Y_{l m}
to real (tesseral) Z_{l m}. Orthogonality is preserved.
"""
function complex_to_real_tensor(Ccx::AbstractArray{<:Number}, ls::Vector{Int}, Lf::Int;
	tol::Real = 1e-12, rephase::Symbol = :global)
	N = length(ls)
	C = Ccx

	# 1) Site-side (bra): apply transformation from complex to real
	for i in 1:N
		S = c2r_sph_harm_matrix(ls[i])  # real = C * complex
		C = nmode_mul(C, conj.(S), i)   # For coefficient tensors, apply C* transformation
	end

	# 2) Final multiplet (ket): apply transformation from complex to real (same as bra side for coefficients)
	Sfinal = c2r_sph_harm_matrix(Lf)
	C = nmode_mul(C, Sfinal, N+1)

	# 3) Phase alignment (rephase). For integer l, a simple analytic global phase
	#    based on the total number of Y_{lm} factors is often sufficient.
	if rephase == :global
		# Use a deterministic global phase instead of data-dependent one:
		# For integer ls, the transformed tensors are guaranteed to be
		# purely real or purely imaginary up to a known factor of i.
		# We compensate that factor here so that the final tensor becomes real.
		# NOTE: This keeps the overall phase convention stable across runs.
		C .*= exp(-1im * (π / 2) * (sum(ls) - Lf))
	elseif rephase == :perM
		# Per-Mf phase alignment: align phase independently for each Mf slice
		allcols = ntuple(_ -> Colon(), N)  # N colons
		for q in axes(C, N+1)
			s = @view C[allcols..., q]
			_, lin = findmax(abs.(s))
			φ = angle(s[lin])
			s .*= exp(-1im * φ)
		end
	end

	# 4) Check imaginary parts and convert to real
	max_imag = maximum(abs.(imag.(C)))
	@assert max_imag ≤ tol "complex_to_real_tensor: imaginary residue too large = $(max_imag)"

	Creal = Array{Float64}(real.(C))  # Always real at this point (type is also Float64)
	return Creal
end

# ---------------- 2') Build complex coefficient tensor with m-indexed axes ----------------

"""
	coeff_tensor_complex_mindexed(ls, Lseq, Lf; T=Float64) -> OffsetArray

	Same coefficients as `coeff_tensor_complex`, but returns an OffsetArray whose axes
	are directly indexed by the magnetic quantum numbers m for each site and Mf for the
	final multiplet. This allows accessing entries as `C[m1, m2, ..., mN, Mf]` with
	`m_k ∈ -l_k: +l_k` and `Mf ∈ -Lf: +Lf`.

	Notes:
	- Computation reuses the same dynamic-programming scheme; only storage layout differs.
	- This is convenient for readability but introduces an OffsetArrays dependency.
	"""
function coeff_tensor_complex_mindexed(
	ls::Vector{Int},
	Lseq::Vector{Int},
	Lf::Int;
	T::Type{<:Real} = Float64,
)
	N = length(ls)
	mr = [mrange(l) for l in ls]                 # per-site m ranges: -l..+l
	dims = [2*l + 1 for l in ls]                 # per-site lengths

	# Allocate a base (0-based) array and wrap with offsets so that indices are actual m values
	Cbase = zeros(T, dims..., 2*Lf + 1)
	site_offsets = Tuple(first.(mr) .- 1)         # shift so that C[m] is valid
	final_offset = (-Lf - 1)                      # for Mf ∈ -Lf..+Lf
	C = OffsetArray(Cbase, site_offsets..., final_offset)

	# N = 1: delta coupling constraint
	if N == 1
		l1 = ls[1]
		for Mf in mrange(Lf)
			# nonzero only when Lf == l1 and Mf == m1
			C[Mf, Mf] = (Lf == l1) ? one(T) : zero(T)
		end
		return C
	end

	# N = 2: single CG coefficient
	if N == 2
		for mi in mr[1], mj in mr[2], Mf in mrange(Lf)
			C[mi, mj, Mf] = convert(T, clebschgordan(ls[1], mi, ls[2], mj, Lf, Mf))
		end
		return C
	end

	# N >= 3: forward dynamic programming over intermediate M's
	Lprev = Lseq[end]
	inds  = CartesianIndices(Tuple(dims))  # Convert Vector to Tuple for multi-dimensional indexing
	for I in inds
		# Map multi-index I to actual magnetic quantum numbers m_k
		ms = ntuple(i -> mr[i][I[i]], N)

		# Stage 1: (l1,l2) -> L12 (sum over M12)
		L12 = Lseq[1]
		acc_prev = Dict{Int, T}()                  # map: M12 -> amplitude
		for M12 in (-L12):L12
			v = convert(T, clebschgordan(ls[1], ms[1], ls[2], ms[2], L12, M12))
			if v != 0
				acc_prev[M12] = v
			end
		end

		# Intermediate stages: (L_{...}, l_{t+1}) -> L_t
		for t in 2:(N-2)
			Lt_1, Lt = Lseq[t-1], Lseq[t]
			mt1 = ms[t+1]
			acc = Dict{Int, T}()
			for (Mprev, aval) in acc_prev
				for Mt in (-Lt):Lt
					v = aval * convert(T, clebschgordan(Lt_1, Mprev, ls[t+1], mt1, Lt, Mt))
					if v != 0
						acc[Mt] = get(acc, Mt, zero(T)) + v
					end
				end
			end
			acc_prev = acc
		end

		# Final stage: (Lprev, lN) -> Lf ; sum over Mprev
		mN = ms[end]
		for Mf in mrange(Lf)
			total_amplitude = zero(T)
			for (Mprev, aval) in acc_prev
				total_amplitude +=
					aval * convert(T, clebschgordan(Lprev, Mprev, ls[end], mN, Lf, Mf))
			end
			C[ms..., Mf] = total_amplitude
		end
	end

	return C
end

# ---------------- 4) Build all bases, grouped by final Lf ----------------


"""
	build_all_complex_bases(ls; normalize=:none, time_reversal::Bool=true, isotropy::Bool=true)
	  -> Dict{Int, Vector{Array{Float64}}}, Dict{Int, Vector{Vector{Int}}}

Build and group all complex-basis coefficient tensors C^(Lf,M) by final Lf.
Return (bases_by_L, paths_by_L):
  - bases_by_L[Lf] = vector of (N+1)-way tensors for each coupling path
  - paths_by_L[Lf] = corresponding Lseq for each tensor

Options:
  normalize:
	:none -> no normalization (pure CG chaining; already orthogonal)
	:fro  -> Frobenius-norm normalization
  isotropy:
	if true (default), keep only Lf == 0 (isotropic scalar). This filtering
	takes precedence over time_reversal.
"""
function build_all_complex_bases(ls::Vector{Int}; normalize::Symbol = :none, isotropy::Bool = false)
	bases_by_L = Dict{Int, Vector{Array{Float64}}}()
	paths_by_L = Dict{Int, Vector{Vector{Int}}}()

	for (Lseq, Lf) in enumerate_paths_left_all(ls)
		if isotropy && Lf != 0
			continue
		end
		Ccx = coeff_tensor_complex(ls, Lseq, Lf)
		if normalize == :fro
			nrm = norm(Ccx)
			if nrm > 0
				Ccx ./= nrm
			end
		elseif normalize != :none
			error("unknown normalize option: $normalize (use :none or :fro)")
		end
		push!(get!(bases_by_L, Lf, []), Array{Float64}(Ccx))
		push!(get!(paths_by_L, Lf, []), Lseq)
	end
	return bases_by_L, paths_by_L
end

"""
	build_all_real_bases(ls; normalize=:none, time_reversal::Bool=true, isotropy::Bool=true)
	  -> Dict{Int, Vector{Array{Float64}}}, Dict{Int, Vector{Vector{Int}}}

Build and group all real (tesseral) coefficient tensors by final Lf.
Same return as build_all_complex_bases, but tensors are already transformed to real
for every site and for the final multiplet axis.

Options:
  normalize:
	:none -> no normalization (orthogonality preserved by unitary transforms)
	:fro  -> Frobenius-norm normalization
"""
function build_all_real_bases(ls::Vector{Int}; normalize::Symbol = :none, isotropy::Bool = false)
	bases_by_L = Dict{Int, Vector{Array{Float64}}}()
	paths_by_L = Dict{Int, Vector{Vector{Int}}}()

	for (Lseq, Lf) in enumerate_paths_left_all(ls)
		if isotropy && Lf != 0
			continue
		end
		Ccx   = coeff_tensor_complex(ls, Lseq, Lf)
		Creal = complex_to_real_tensor(Ccx, ls, Lf)
		if normalize == :fro
			nrm = norm(Creal)
			if nrm > 0
				Creal ./= nrm
			end
		elseif normalize != :none
			error("unknown normalize option: $normalize (use :none or :fro)")
		end
		push!(get!(bases_by_L, Lf, []), Creal)
		push!(get!(paths_by_L, Lf, []), Lseq)
	end
	return bases_by_L, paths_by_L
end

# ---------------- Exports ----------------

export mrange,
	enumerate_paths_left_all,
	nmode_mul,
	coeff_tensor_complex,
	coeff_tensor_complex_mindexed,
	complex_to_real_tensor,
	build_all_complex_bases,
	build_all_real_bases

end # module
