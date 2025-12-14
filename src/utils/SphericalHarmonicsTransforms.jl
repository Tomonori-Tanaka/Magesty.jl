module SphericalHarmonicsTransforms

"""
    c2r_sph_harm_matrix(l::Integer) -> Matrix{ComplexF64}

Compute the transformation matrix from complex spherical harmonics to real spherical harmonics.
"""
function c2r_sph_harm_matrix(l::Integer)::Matrix{ComplexF64}
    l ≥ 0 || throw(ArgumentError("Angular momentum l must be non-negative"))
    
    U = zeros(ComplexF64, 2l+1, 2l+1)
    idx(m) = m + l + 1

    # m = 0
    U[idx(0), idx(0)] = 1

    # m > 0 blocks (pair of real columns for ±m)
    for m in 1:l
        jp, jm = idx(+m), idx(-m)
        # Row for complex +m
        U[idx(+m), jp] = (-1)^m / sqrt(2)
        U[idx(+m), jm] = im * (-1)^m / sqrt(2)
        # Row for complex -m
        U[idx(-m), jp] = 1 / sqrt(2)
        U[idx(-m), jm] = -im / sqrt(2)
    end
    return U
end


"""
    r2c_sph_harm_matrix(l::Integer) -> Matrix{ComplexF64}

Compute the transformation matrix from real spherical harmonics to complex spherical harmonics.
"""
function r2c_sph_harm_matrix(l::Integer)::Matrix{ComplexF64}
    l ≥ 0 || throw(ArgumentError("Angular momentum l must be non-negative"))
    
    U = zeros(ComplexF64, 2l + 1, 2l + 1)
    idx(m) = m + l + 1

    # m = 0
    U[idx(0), idx(0)] = 1

    # m > 0 blocks (pair of real columns for ±m)
    for m in 1:l
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