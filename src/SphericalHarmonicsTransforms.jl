module SphericalHarmonicsTransforms

export r2c_sph_harm_matrix, c2r_sph_harm_matrix

"""
    r2c_sph_harm_matrix(l::Integer) -> Matrix{ComplexF64}

Compute the transformation matrix from tesseral harmonics to complex spherical harmonics.

This is the inverse of [`c2r_sph_harm_matrix`](@ref), obtained as its adjoint.
The adjoint equals the inverse only because `c2r_sph_harm_matrix` is unitary
(an `L2`-orthonormal change of basis between the complex and real harmonics);
a future change to the normalization that breaks that unitarity would break
this identity. The round-trip `r2c * c2r ≈ I` is asserted in the tests.
"""
function r2c_sph_harm_matrix(l::Integer)::Matrix{ComplexF64}
    # Adjoint = inverse because c2r is unitary (see the docstring note).
    U = c2r_sph_harm_matrix(l)
    U = U'
    return U
end


"""
    c2r_sph_harm_matrix(l::Integer) -> Matrix{ComplexF64}

Compute the transformation matrix from complex spherical harmonics to tesseral harmonics.
"""
function c2r_sph_harm_matrix(l::Integer)::Matrix{ComplexF64}
    l ≥ 0 || throw(ArgumentError("Angular momentum l must be non-negative"))
    
    U = zeros(ComplexF64, 2l + 1, 2l + 1)
    idx(m) = m + l + 1

    # m = 0
    U[idx(0), idx(0)] = 1

    # m > 0 blocks (pair of real columns for ±m)
    for m = 1:l
        ip, im_idx = idx(+m), idx(-m)
        # Row for complex +m
        U[idx(+m), ip] = (-1)^m / sqrt(2)
        U[idx(+m), im_idx] = 1 / sqrt(2)
        # Row for complex -m
        U[idx(-m), ip] = -im * (-1)^m / sqrt(2)
        U[idx(-m), im_idx] = im / sqrt(2)
    end
    return U
end

end # module SphericalHarmonicsTransforms