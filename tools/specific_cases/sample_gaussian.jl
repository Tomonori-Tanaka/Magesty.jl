using Random
using LinearAlgebra

"""
    sample_spherical_gaussian(mu::Vector{Float64}, sigma::Float64)

Generates a random point on the sphere following a Gaussian-like distribution.
The distribution is centered at direction `mu` (a 3D unit vector) and has an angular dispersion `sigma` (in radians).
The probability density for the great circle distance θ is proportional to
    f(θ) ∝ exp( -θ²/(2σ²) ) * sinθ.
"""
function sample_spherical_gaussian(mu::Vector{<:Real}, sigma::Real)::Vector{Float64}
    # Gaussian-weighted density for theta (normalization constant is not needed)
    f(theta) = sin(theta) * exp(-theta^2 / (2*sigma^2))
    
    # Determine the maximum of f(theta) on [0, π] for rejection sampling
    N_grid = 1000
    theta_vals = range(0, stop=pi, length=N_grid)
    M_val = maximum(f.(theta_vals))
    
    # Rejection sampling for theta in [0, π]
    theta = 0.0
    while true
        theta_candidate = rand() * pi  # uniform sample in [0, π]
        if rand() * M_val < f(theta_candidate)
            theta = theta_candidate
            break
        end
    end
    
    # Sample phi uniformly from [0, 2π)
    phi = 2*pi * rand()
    
    # Convert spherical coordinates (centered at the north pole) to Cartesian coordinates
    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)
    p = [x, y, z]
    
    # Rotate the sample point from the north pole ([0,0,1]) to the specified direction mu.
    north = [0.0, 0.0, 1.0]
    if norm(mu - north) < 1e-10
        # Already aligned with the north pole; return the sample directly.
        return p
    elseif norm(mu + north) < 1e-10
        # If mu is the opposite of the north pole, apply a 180° rotation around the x-axis.
        R = [1.0  0.0  0.0;
             0.0 -1.0  0.0;
             0.0  0.0 -1.0]
        return R * p
    else
        # Determine the rotation axis and angle.
        axis = cross(north, mu)
        axis = axis / norm(axis)
        angle = acos(clamp(dot(north, mu), -1.0, 1.0))
        # Construct the rotation matrix using Rodrigues' rotation formula.
        K = [  0.0      -axis[3]   axis[2];
             axis[3]    0.0      -axis[1];
            -axis[2]   axis[1]    0.0    ]
        R = I + sin(angle)*K + (1-cos(angle))*(K*K)
        return R * p
    end
end
