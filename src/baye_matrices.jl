"""
    Baye's exact D and T matrices for Lagrange-Legendre basis.

These matrices act on expansion coefficients c_j, NOT on function values φ(x_j).

For x-regularized basis: f_j(r) = α_j * r * L_j(r/R)
Wave function: φ(r) = Σ_j c_j * f_j(r)
At mesh points: φ(r_j) = c_j * α_j * r_j
So: c_j = φ(r_j) / (α_j * r_j)

Reference: D. Baye, Physics Reports 565 (2015) 1-107, Section 3.4.5
           Equations (3.122)-(3.126)
"""

"""
    baye_D_matrix(mesh::LagrangeMesh) -> Matrix{Float64}

Compute Baye's exact D matrix (first derivative d/dx on unit interval).

D_{i≠j} = (-1)^{i-j} * √[x_i(1-x_j)/(x_j(1-x_i))] / (x_i - x_j)
D_{ii} = 1 / [2*x_i*(1-x_i)]

The returned matrix is scaled to physical coordinates: D → D/R

# Arguments
- `mesh::LagrangeMesh`: The Lagrange mesh

# Returns
- `Matrix{Float64}`: N×N derivative matrix scaled by 1/R
"""
function baye_D_matrix(mesh::LagrangeMesh)
    N = mesh.N
    x = mesh.x
    R = mesh.R

    D = zeros(Float64, N, N)

    for i in 1:N
        for j in 1:N
            if i != j
                # Off-diagonal: D_{ij} = (-1)^{i-j} * sqrt[x_i(1-x_j)/(x_j(1-x_i))] / (x_i - x_j)
                sign = (-1)^(i - j)
                num = x[i] * (1 - x[j])
                den = x[j] * (1 - x[i])
                D[i, j] = sign * sqrt(num / den) / (x[i] - x[j])
            else
                # Diagonal: D_{ii} = 1 / [2*x_i*(1-x_i)]
                D[i, i] = 1.0 / (2 * x[i] * (1 - x[i]))
            end
        end
    end

    # Scale to physical interval (0, R)
    return D / R
end

"""
    baye_T_matrix(mesh::LagrangeMesh) -> Matrix{Float64}

Compute Baye's exact T matrix (kinetic energy: -d²/dx² on unit interval).

For Lagrange-Legendre on (0,1):

T_{i≠j} = (-1)^{i-j} * (x_i + x_j - 2*x_i²) / [x_j*(x_j-x_i)²]
          * √[x_j*(1-x_j) / (x_i*(1-x_i)³)]

T_{ii} = [N(N+1)*x_i*(1-x_i) - 3*x_i + 1] / [3*x_i²*(1-x_i)²]

The returned matrix is scaled to physical coordinates: T → T/R²

Note: This is the NEGATIVE of the second derivative: T = -d²/dx²

# Arguments
- `mesh::LagrangeMesh`: The Lagrange mesh

# Returns
- `Matrix{Float64}`: N×N kinetic energy matrix scaled by 1/R²
"""
function baye_T_matrix(mesh::LagrangeMesh)
    N = mesh.N
    x = mesh.x
    R = mesh.R

    T = zeros(Float64, N, N)

    for i in 1:N
        xi = x[i]
        for j in 1:N
            xj = x[j]
            if i != j
                # Off-diagonal element
                sign = (-1)^(i - j)

                num = xi + xj - 2 * xi^2
                den = xj * (xj - xi)^2

                sqrt_num = xj * (1 - xj)
                sqrt_den = xi * (1 - xi)^3

                T[i, j] = sign * (num / den) * sqrt(sqrt_num / sqrt_den)
            else
                # Diagonal element
                # T_{ii} = [N(N+1)*x_i*(1-x_i) - 3*x_i + 1] / [3*x_i²*(1-x_i)²]
                num = N * (N + 1) * xi * (1 - xi) - 3 * xi + 1
                den = 3 * xi^2 * (1 - xi)^2
                T[i, i] = num / den
            end
        end
    end

    # Scale to physical interval (0, R)
    return T / R^2
end

"""
    kinetic_matrix(mesh::LagrangeMesh) -> Matrix{Float64}

Compute the kinetic energy matrix in the Lagrange basis.
This is essentially the T matrix for use in the Hamiltonian.
"""
function kinetic_matrix(mesh::LagrangeMesh)
    return baye_T_matrix(mesh)
end
