"""
    SLAM.jl - Scattering Lagrange Asymptotic Matching

A Julia package for solving nuclear scattering problems using Lagrange mesh methods
with exact differential matrices operating on expansion coefficients.

Key features:
- Lagrange-Legendre mesh on (0, R) using Gauss-Legendre quadrature
- Analytical formulas for first and second derivative matrices (D and T)
- x-regularized basis functions: f_j(r) = α_j * r * L_j(r/R)
- Coefficient-space formulation for proper matrix action
- Outgoing wave boundary condition using Coulomb-Hankel functions
- COUL90 Fortran library integration for accurate Coulomb functions
- Woods-Saxon optical model potentials with spin-orbit coupling

## References
- D. Baye, Physics Reports 565 (2015) 1-107
- A.R. Barnett, CPC 21 (1981) 297-314 (COUL90)

## Usage
```julia
using SLAM

# Define optical potential
pot = OpticalPotential(
    V_v=50.0, r_v=1.25, a_v=0.65,    # Real volume
    W_v=10.0, r_wv=1.25, a_wv=0.65,  # Imaginary volume
    r_c=1.25,                         # Coulomb radius
    Z_proj=1.0, A_proj=1.0,          # Proton
    Z_targ=20.0, A_targ=40.0         # Ca-40
)

# Create scattering problem
prob = ScatteringProblem(pot, 30.0, 0)  # E=30 MeV, l=0

# Solve
result = solve_scattering(prob)
println("S-matrix: ", result.S_matrix)
println("Cross section: ", result.cross_section, " mb")
```
"""
module SLAM

using LinearAlgebra
using FastGaussQuadrature
using SpecialFunctions
using Printf

# Include submodules
include("mesh.jl")
include("baye_matrices.jl")
include("coulomb.jl")
include("potentials.jl")
include("scattering.jl")
include("numerov.jl")

# Module initialization
function __init__()
    __init_coulomb__()
end

# Export mesh functions
export LagrangeMesh, init_legendre_mesh
export lagrange_function, lagrange_derivative
export basis_function_at_R, basis_derivative_at_R
export basis_function_at_R_analytical, basis_derivative_at_R_analytical
export fhat_at_boundary, dfhat_dx_at_boundary
# Transformed mesh exports
export MeshTransform, TRANSFORM_IDENTITY, TRANSFORM_POWER, TRANSFORM_FRESCO
export TransformedLagrangeMesh, init_transformed_mesh
export mesh_point_distribution

# Export Baye matrix functions
export baye_D_matrix, baye_T_matrix, kinetic_matrix
export numerical_D_matrix, numerical_T_matrix

# Export Coulomb functions
export coulomb_eta, coulomb_k
export coulomb_F, coulomb_G, coulomb_FG
export coulomb_H_plus, coulomb_H_plus_derivative
export coulomb_gamma_s, coulomb_sigma

# Export potential functions
export OpticalPotential, reduced_mass
export woods_saxon, woods_saxon_derivative
export evaluate_potential, evaluate_coulomb
export evaluate_short_range, evaluate_total_potential

# Export scattering functions
export ScatteringProblem, ScatteringResult
export SourceTermMethod, LAGRANGE_GAUSS, FINE_GRID
export solve_scattering, SMatrix, cross_section
export solve_all_partial_waves
export elastic_differential_cross_section
export legendre_P
export compute_source_term_fine

# Export wave function functions
export WaveFunctionResult
export solve_wavefunction, evaluate_wavefunction

# Export Numerov functions
export NumerovResult
export solve_numerov, compare_numerov_lagrange

end # module SLAM
