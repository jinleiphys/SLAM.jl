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
