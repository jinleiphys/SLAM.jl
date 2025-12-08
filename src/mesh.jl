"""
    Lagrange-Legendre mesh implementation for scattering calculations.

Uses Gauss-Legendre quadrature points on (0,1) scaled to (0,R).

The x-regularized basis functions are:
f_j(r) = (r/r_j) * L_j(r) / √λ_j

where L_j(r) is the standard Lagrange polynomial satisfying L_j(r_k) = δ_{jk}

At mesh points: f_j(r_j) = 1/√λ_j
So: c_j = φ(r_j) * √λ_j

Reference: D. Baye, Physics Reports 565 (2015) 1-107, Section 3.4
"""

"""
    LagrangeMesh

Structure holding Lagrange-Legendre mesh data.

# Fields
- `N::Int`: Number of mesh points
- `R::Float64`: Outer boundary radius
- `x::Vector{Float64}`: Mesh points in (0,1)
- `r::Vector{Float64}`: Physical mesh points in (0,R)
- `w::Vector{Float64}`: Gauss-Legendre weights on (0,1)
- `λ::Vector{Float64}`: Scaled weights (λ_j = w_j/2)
- `α::Vector{Float64}`: Normalization factors (α_j = 1/√λ_j)
"""
struct LagrangeMesh
    N::Int
    R::Float64
    x::Vector{Float64}   # Mesh points in (0,1)
    r::Vector{Float64}   # Physical mesh points r = R*x
    w::Vector{Float64}   # Gauss weights on (0,1)
    λ::Vector{Float64}   # Scaled weights λ_j = w_j/2
    α::Vector{Float64}   # Normalization α_j = 1/√λ_j
end

"""
    init_legendre_mesh(N::Int, R::Float64) -> LagrangeMesh

Initialize a Lagrange-Legendre mesh with N points on the interval (0, R).
"""
function init_legendre_mesh(N::Int, R::Float64)
    # Get Gauss-Legendre nodes and weights on (-1, 1)
    t, w_gl = gausslegendre(N)

    # Transform to (0, 1)
    x = @. (t + 1) / 2
    w = @. w_gl / 2

    # Scale to physical coordinates
    r = @. R * x

    # Compute λ and α
    λ = copy(w)
    α = @. 1.0 / sqrt(λ)

    return LagrangeMesh(N, R, x, r, w, λ, α)
end

"""
    lagrange_polynomial(mesh::LagrangeMesh, j::Int, r_eval::Float64) -> Float64

Evaluate the j-th Lagrange polynomial at r_eval (in physical coordinates).

L_j(r) = ∏_{k≠j} (r - r_k) / (r_j - r_k)

Note: L_j(r_k) = δ_{jk}
"""
function lagrange_polynomial(mesh::LagrangeMesh, j::Int, r_eval::Float64)
    result = 1.0
    for k in 1:mesh.N
        if k != j
            result *= (r_eval - mesh.r[k]) / (mesh.r[j] - mesh.r[k])
        end
    end
    return result
end

"""
    lagrange_polynomial_derivative(mesh::LagrangeMesh, j::Int, r_eval::Float64) -> Float64

Evaluate the derivative of the j-th Lagrange polynomial at r_eval (physical coords).
"""
function lagrange_polynomial_derivative(mesh::LagrangeMesh, j::Int, r_eval::Float64)
    result = 0.0
    for m in 1:mesh.N
        if m != j
            term = 1.0 / (mesh.r[j] - mesh.r[m])
            for k in 1:mesh.N
                if k != j && k != m
                    term *= (r_eval - mesh.r[k]) / (mesh.r[j] - mesh.r[k])
                end
            end
            result += term
        end
    end
    return result
end

"""
    basis_function_at_R(mesh::LagrangeMesh, j::Int) -> Float64

Evaluate the x-regularized basis function f_j(r) at r = R.

f_j(r) = (r/r_j) * L_j(r) / √λ_j

At r = R:
f_j(R) = (R/r_j) * L_j(R) / √λ_j
"""
function basis_function_at_R(mesh::LagrangeMesh, j::Int)
    L_j_at_R = lagrange_polynomial(mesh, j, mesh.R)
    return (mesh.R / mesh.r[j]) * L_j_at_R / sqrt(mesh.λ[j])
end

"""
    basis_derivative_at_R(mesh::LagrangeMesh, j::Int) -> Float64

Evaluate the derivative of the x-regularized basis function at r = R.

f_j(r) = (r/r_j) * L_j(r) / √λ_j

f_j'(r) = (1/r_j) * L_j(r) / √λ_j + (r/r_j) * L_j'(r) / √λ_j
        = [L_j(r)/r_j + (r/r_j)*L_j'(r)] / √λ_j

At r = R:
f_j'(R) = [L_j(R)/r_j + (R/r_j)*L_j'(R)] / √λ_j
"""
function basis_derivative_at_R(mesh::LagrangeMesh, j::Int)
    L_j_at_R = lagrange_polynomial(mesh, j, mesh.R)
    dL_j_at_R = lagrange_polynomial_derivative(mesh, j, mesh.R)
    return (L_j_at_R / mesh.r[j] + (mesh.R / mesh.r[j]) * dL_j_at_R) / sqrt(mesh.λ[j])
end

"""
    fhat_at_boundary(mesh::LagrangeMesh, j::Int) -> Float64

Evaluate the x-regularized basis function f̂_j(x) at x = 1 (r = R) using analytical formula.

From DBMM paper (Eq. 259), using P_N(1) = 1:
    f̂_j(1) = (-1)^{N-j} / √[x_j(1-x_j)]

This is the value of f̂_j defined in Eq. 195:
    f̂_j(x) = (-1)^{N-j} √[(1-x_j)/x_j] · x·P_N(2x-1)/(x - x_j)

These satisfy f̂_j(x_i) = δ_{ij}/√λ_j
"""
function fhat_at_boundary(mesh::LagrangeMesh, j::Int)
    x_j = mesh.x[j]
    sign = (-1)^(mesh.N - j)
    return sign / sqrt(x_j * (1 - x_j))
end

"""
    dfhat_dx_at_boundary(mesh::LagrangeMesh, j::Int) -> Float64

Evaluate df̂_j/dx at x = 1 using analytical formula.

From DBMM paper (Eq. 264), using P_N'(1) = N(N+1)/2:
    df̂_j/dx|_{x=1} = (-1)^{N-j}/√[x_j(1-x_j)] * [N(N+1) - x_j/(1-x_j)]
"""
function dfhat_dx_at_boundary(mesh::LagrangeMesh, j::Int)
    x_j = mesh.x[j]
    N = mesh.N
    sign = (-1)^(N - j)

    prefactor = sign / sqrt(x_j * (1 - x_j))
    bracket = N * (N + 1) - x_j / (1 - x_j)

    return prefactor * bracket
end

"""
    basis_function_at_R_analytical(mesh::LagrangeMesh, j::Int) -> Float64

For backward compatibility - returns f̂_j(1).
"""
function basis_function_at_R_analytical(mesh::LagrangeMesh, j::Int)
    return fhat_at_boundary(mesh, j)
end

"""
    basis_derivative_at_R_analytical(mesh::LagrangeMesh, j::Int) -> Float64

For backward compatibility - returns df̂_j/dx at x=1.
"""
function basis_derivative_at_R_analytical(mesh::LagrangeMesh, j::Int)
    return dfhat_dx_at_boundary(mesh, j)
end

# Keep old names for compatibility
lagrange_function(mesh, j, x_eval) = lagrange_polynomial(mesh, j, x_eval * mesh.R) * sqrt(mesh.λ[j])
lagrange_derivative(mesh, j, x_eval) = lagrange_polynomial_derivative(mesh, j, x_eval * mesh.R) * mesh.R * sqrt(mesh.λ[j])

#==============================================================================
# Coordinate-transformed Lagrange mesh
#
# Uses Gauss-Legendre points on (0,1) but applies a nonlinear transformation
# r = R * g(x) to concentrate mesh points in the small-r region where
# the potential and source term are significant.
#
# Available transformations:
# - :identity  -> g(x) = x        (standard Legendre mesh)
# - :power     -> g(x) = x^p      (p > 1 concentrates points near r=0)
# - :fresco    -> g(x) = (3x²+1)x/4  (similar to FRESCO code)
#
# Reference: Thompson & Nunes, "Nuclear Reactions for Astrophysics", Eq. 133
==============================================================================#

"""
    MeshTransform

Enum for mesh coordinate transformation types.
"""
@enum MeshTransform begin
    TRANSFORM_IDENTITY = 1   # g(x) = x (standard Legendre)
    TRANSFORM_POWER = 2      # g(x) = x^p
    TRANSFORM_FRESCO = 3     # g(x) = (3x²+1)x/4
end

"""
    TransformedLagrangeMesh

Structure holding Lagrange mesh data with coordinate transformation.

The transformation r = R * g(x) maps Legendre points x ∈ (0,1) to physical
coordinates r ∈ (0,R), concentrating points in the small-r region.

# Fields
- `N::Int`: Number of mesh points
- `R::Float64`: Outer boundary radius
- `x::Vector{Float64}`: Original Legendre points in (0,1)
- `r::Vector{Float64}`: Physical mesh points after transformation
- `w::Vector{Float64}`: Gauss-Legendre weights on (0,1)
- `w_phys::Vector{Float64}`: Integration weights in physical coordinates
- `λ::Vector{Float64}`: Weights for basis function normalization
- `α::Vector{Float64}`: Normalization factors (α_j = 1/√λ_j)
- `g::Vector{Float64}`: Transformation function values g(x_j)
- `dg::Vector{Float64}`: Transformation derivative g'(x_j)
- `transform::MeshTransform`: Type of transformation used
- `power::Float64`: Power parameter (only for TRANSFORM_POWER)
"""
struct TransformedLagrangeMesh
    N::Int
    R::Float64
    x::Vector{Float64}      # Original Legendre points (0,1)
    r::Vector{Float64}      # Physical coordinates
    w::Vector{Float64}      # Gauss weights on (0,1)
    w_phys::Vector{Float64} # Physical integration weights
    λ::Vector{Float64}      # Normalization weights
    α::Vector{Float64}      # Normalization factors
    g::Vector{Float64}      # g(x) values
    dg::Vector{Float64}     # g'(x) values
    transform::MeshTransform
    power::Float64
end

"""
    init_transformed_mesh(N::Int, R::Float64; transform=TRANSFORM_IDENTITY, power=2.0)

Initialize a Lagrange mesh with coordinate transformation.

The transformation r = R * g(x) concentrates mesh points near r = 0.

# Arguments
- `N::Int`: Number of mesh points
- `R::Float64`: Outer boundary radius
- `transform::MeshTransform`: Type of transformation (default: TRANSFORM_IDENTITY)
- `power::Float64`: Power for TRANSFORM_POWER (default: 2.0)

# Transformations
- `TRANSFORM_IDENTITY`: g(x) = x (standard Legendre mesh)
- `TRANSFORM_POWER`: g(x) = x^p where p = power parameter
- `TRANSFORM_FRESCO`: g(x) = (3x² + 1)x / 4 (concentrates ~22% of points in r < R/30)

# Returns
- `TransformedLagrangeMesh`: Mesh structure with transformed coordinates

# Example
```julia
# Standard Legendre mesh
mesh = init_transformed_mesh(40, 30.0)

# Power transformation with p=2 (28% of points in r < R/30)
mesh = init_transformed_mesh(40, 30.0, transform=TRANSFORM_POWER, power=2.0)

# FRESCO-style transformation
mesh = init_transformed_mesh(40, 30.0, transform=TRANSFORM_FRESCO)
```
"""
function init_transformed_mesh(N::Int, R::Float64;
                                transform::MeshTransform=TRANSFORM_IDENTITY,
                                power::Float64=2.0)
    # Get Gauss-Legendre nodes and weights on (-1, 1)
    t, w_gl = gausslegendre(N)

    # Transform to (0, 1)
    x = @. (t + 1) / 2
    w = @. w_gl / 2

    # Apply coordinate transformation r = R * g(x)
    if transform == TRANSFORM_IDENTITY
        g = copy(x)
        dg = ones(N)
    elseif transform == TRANSFORM_POWER
        g = x.^power
        dg = power .* x.^(power - 1)
    elseif transform == TRANSFORM_FRESCO
        g = @. (3*x^2 + 1) * x / 4
        dg = @. (9*x^2 + 1) / 4
    else
        error("Unknown transform: $transform")
    end

    r = R .* g

    # Physical integration weights: ∫f(r)dr = ∫f(R*g(x)) * R*g'(x) dx
    w_phys = w .* R .* dg

    # For basis function normalization, use the original Legendre weights
    λ = copy(w)
    α = @. 1.0 / sqrt(λ)

    return TransformedLagrangeMesh(N, R, x, r, w, w_phys, λ, α, g, dg, transform, power)
end

# Make TransformedLagrangeMesh compatible with existing functions
# by implementing the same interface as LagrangeMesh

function lagrange_polynomial(mesh::TransformedLagrangeMesh, j::Int, r_eval::Float64)
    result = 1.0
    for k in 1:mesh.N
        if k != j
            result *= (r_eval - mesh.r[k]) / (mesh.r[j] - mesh.r[k])
        end
    end
    return result
end

function lagrange_polynomial_derivative(mesh::TransformedLagrangeMesh, j::Int, r_eval::Float64)
    result = 0.0
    for m in 1:mesh.N
        if m != j
            term = 1.0 / (mesh.r[j] - mesh.r[m])
            for k in 1:mesh.N
                if k != j && k != m
                    term *= (r_eval - mesh.r[k]) / (mesh.r[j] - mesh.r[k])
                end
            end
            result += term
        end
    end
    return result
end

function basis_function_at_R(mesh::TransformedLagrangeMesh, j::Int)
    L_j_at_R = lagrange_polynomial(mesh, j, mesh.R)
    return (mesh.R / mesh.r[j]) * L_j_at_R / sqrt(mesh.λ[j])
end

function basis_derivative_at_R(mesh::TransformedLagrangeMesh, j::Int)
    L_j_at_R = lagrange_polynomial(mesh, j, mesh.R)
    dL_j_at_R = lagrange_polynomial_derivative(mesh, j, mesh.R)
    return (L_j_at_R / mesh.r[j] + (mesh.R / mesh.r[j]) * dL_j_at_R) / sqrt(mesh.λ[j])
end

"""
    fhat_at_boundary(mesh::TransformedLagrangeMesh, j::Int) -> Float64

Evaluate the x-regularized basis function f̂_j(x) at x = 1 (r = R) using analytical formula.
Same as for LagrangeMesh since basis functions are defined in x-space.
"""
function fhat_at_boundary(mesh::TransformedLagrangeMesh, j::Int)
    x_j = mesh.x[j]
    sign = (-1)^(mesh.N - j)
    return sign / sqrt(x_j * (1 - x_j))
end

"""
    dfhat_dx_at_boundary(mesh::TransformedLagrangeMesh, j::Int) -> Float64

Evaluate df̂_j/dx at x = 1 using analytical formula.
Same as for LagrangeMesh since basis functions are defined in x-space.
"""
function dfhat_dx_at_boundary(mesh::TransformedLagrangeMesh, j::Int)
    x_j = mesh.x[j]
    N = mesh.N
    sign = (-1)^(N - j)
    prefactor = sign / sqrt(x_j * (1 - x_j))
    bracket = N * (N + 1) - x_j / (1 - x_j)
    return prefactor * bracket
end

"""
    mesh_point_distribution(mesh::Union{LagrangeMesh, TransformedLagrangeMesh};
                            regions=[(0,1), (1,3), (3,6), (6,10), (10,15), (15,Inf)])

Print statistics about mesh point distribution in different radial regions.
"""
function mesh_point_distribution(mesh::Union{LagrangeMesh, TransformedLagrangeMesh};
                                  regions=nothing)
    R = mesh.R
    if regions === nothing
        # Default regions based on R
        regions = [(0, R/30), (R/30, R/10), (R/10, R/5), (R/5, R/3), (R/3, R/2), (R/2, R)]
    end

    println("Mesh point distribution (N=$(mesh.N), R=$R):")
    println("-"^50)
    for (r_min, r_max) in regions
        r_max_eff = min(r_max, R)
        n = count(r -> r_min <= r < r_max_eff, mesh.r)
        @printf("  [%6.2f, %6.2f) fm: %3d points (%5.1f%%)\n",
                r_min, r_max_eff, n, 100*n/mesh.N)
    end
    println("-"^50)
end
