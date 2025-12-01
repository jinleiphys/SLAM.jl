"""
    Main scattering solver using Baye's method with Lagrange-Legendre basis.

Implements Method 5: Baye's exact D and T matrices acting on expansion coefficients.

Key insight: Baye's D and T matrices act on expansion coefficients c_j,
NOT on function values φ(x_j)!

For x-regularized basis: f_j(r) = (r/r_j) * L_j(r) / √λ_j
Wave function: φ(r) = Σ_j c_j * f_j(r)
At mesh point r_j: f_j(r_j) = 1/√λ_j
So: c_j = φ(r_j) * √λ_j

Reference: D. Baye, Physics Reports 565 (2015) 1-107, Section 3.4.5
"""

# Physical constants
const HBARC = 197.327  # MeV·fm

"""
    ScatteringProblem

Structure holding all parameters for a scattering calculation.
"""
struct ScatteringProblem
    pot::OpticalPotential
    E_cm::Float64
    l::Int
    j::Float64
    N::Int
    R::Float64
end

"""
    ScatteringProblem(pot, E_cm, l; j=nothing, N=60, R=20.0)

Construct a scattering problem with automatic j determination.
"""
function ScatteringProblem(pot::OpticalPotential, E_cm::Float64, l::Int;
                           j::Union{Float64,Nothing}=nothing, N::Int=60, R::Float64=20.0)
    if j === nothing
        j = Float64(l) + 0.5
    end
    return ScatteringProblem(pot, E_cm, l, j, N, R)
end

"""
    ScatteringResult

Structure holding the results of a scattering calculation.
"""
struct ScatteringResult
    S_matrix::ComplexF64
    phase_shift::ComplexF64
    f_amplitude::ComplexF64
    cross_section::Float64
    converged::Bool
end

"""
    solve_scattering(prob::ScatteringProblem) -> ScatteringResult

Solve the scattering problem using Method 5 (Baye's coefficient-space formulation).

The radial Schrödinger equation:
[-ℏ²/(2μ) d²/dr² + ℏ²l(l+1)/(2μr²) + V(r)]φ = Eφ

Rewriting with U = 2μV/ℏ² and k² = 2μE/ℏ²:
[-d²/dr² + l(l+1)/r² + U(r) - k²]φ = 0

For inhomogeneous equation (scattering):
[-d²/dr² + l(l+1)/r² + U(r) - k²]φ = U(r) * F_l(kr)

In matrix form (Method 5 with T = -d²/dr²):
[-T + (k² - l(l+1)/r² - U)] c = U * F_l * √λ
"""
function solve_scattering(prob::ScatteringProblem)
    pot = prob.pot
    E_cm = prob.E_cm
    l = prob.l
    j = prob.j
    N = prob.N
    R = prob.R

    # Compute reduced mass and wave number
    μ = reduced_mass(pot)
    k = coulomb_k(μ, E_cm)
    η = coulomb_eta(pot.Z_proj, pot.Z_targ, μ, E_cm)

    # Conversion factor: U = 2μV/ℏ² (converts V in MeV to U in fm⁻²)
    coeff_kin = HBARC^2 / (2.0 * μ)  # ℏ²/(2μ) in MeV·fm²

    # Initialize mesh
    mesh = init_legendre_mesh(N, R)

    # Compute Baye's T matrix
    T = baye_T_matrix(mesh)

    # Compute boundary condition parameter γ_s = (d/dr)[H⁺(kr)] / H⁺(kR)
    ρ_R = k * R
    γ_s = coulomb_gamma_s(l, η, ρ_R, k)

    # Build the matrix equation M * c = b
    M = zeros(ComplexF64, N, N)
    b = zeros(ComplexF64, N)

    # Interior points (i = 1 to N-1)
    for i in 1:N-1
        r_i = mesh.r[i]
        ρ_i = k * r_i

        # Get potential in MeV and convert to U = 2μV/ℏ² (fm⁻²)
        V_i = evaluate_short_range(pot, r_i, l, j)
        U_i = V_i / coeff_kin

        # Evaluate Coulomb function F_l at ρ_i
        F_l_i = coulomb_F(l, η, ρ_i)

        # Centrifugal term
        cent = Float64(l * (l + 1)) / r_i^2

        # Build matrix row
        # M = -T + (k² - l(l+1)/r² - U) on diagonal
        for jj in 1:N
            M[i, jj] = -T[i, jj]
            if i == jj
                M[i, jj] += k^2 - cent - U_i
            end
        end

        # Source term: U * F_l * √λ
        b[i] = U_i * F_l_i * sqrt(mesh.λ[i])
    end

    # Boundary condition at r = R (row N)
    # φ'(R) - γ_s * φ(R) = 0
    # In basis: Σ_j c_j * [f_j'(R) - γ_s * f_j(R)] = 0
    for jj in 1:N
        f_j_R = basis_function_at_R(mesh, jj)
        f_j_prime_R = basis_derivative_at_R(mesh, jj)
        M[N, jj] = f_j_prime_R - γ_s * f_j_R
    end
    b[N] = 0.0

    # Solve the linear system
    c = M \ b

    # Extract scattering amplitude
    # φ(R) = Σ_j c_j * f_j(R)
    φ_R = sum(c[jj] * basis_function_at_R(mesh, jj) for jj in 1:N)

    # H⁺_l(kR) = G_l + i*F_l
    H_plus_R = coulomb_H_plus(l, η, ρ_R)

    # Scattering amplitude: f_l = φ(R) / [k * H⁺_l(kR)]
    f_l = φ_R / (k * H_plus_R)

    # S-matrix: S_l = 1 + 2ik*f_l
    S_l = 1.0 + 2.0im * k * f_l

    # Phase shift: S_l = exp(2iδ_l)
    δ_l = -0.5im * log(S_l)

    # Partial wave cross section
    # σ_reaction = (π/k²) * (2J+1)/(2S+1) * (1 - |S_l|²) * 10 (fm² to mb)
    s_spin = 0.5  # spin-1/2 projectile
    σ_l = π / k^2 / (2*s_spin + 1) * (2*j + 1) * (1.0 - abs(S_l)^2) * 10.0

    return ScatteringResult(S_l, δ_l, f_l, σ_l, true)
end

"""
    SMatrix(prob::ScatteringProblem) -> ComplexF64

Convenience function to get only the S-matrix element.
"""
function SMatrix(prob::ScatteringProblem)
    result = solve_scattering(prob)
    return result.S_matrix
end

"""
    cross_section(prob::ScatteringProblem) -> Float64

Convenience function to get only the cross section.
"""
function cross_section(prob::ScatteringProblem)
    result = solve_scattering(prob)
    return result.cross_section
end

"""
    solve_all_partial_waves(pot::OpticalPotential, E_cm::Float64;
                            l_max::Int=20, N::Int=60, R::Float64=20.0)

Solve scattering for all partial waves up to l_max.
"""
function solve_all_partial_waves(pot::OpticalPotential, E_cm::Float64;
                                  l_max::Int=20, N::Int=60, R::Float64=20.0)
    results = ScatteringResult[]
    σ_total = 0.0

    for l in 0:l_max
        if l == 0
            j_values = [0.5]
        else
            j_values = [Float64(l) - 0.5, Float64(l) + 0.5]
        end

        for j in j_values
            prob = ScatteringProblem(pot, E_cm, l; j=j, N=N, R=R)
            result = solve_scattering(prob)
            push!(results, result)
            σ_total += result.cross_section
        end
    end

    return results, σ_total
end

"""
    elastic_differential_cross_section(pot::OpticalPotential, E_cm::Float64, θ::Float64;
                                        l_max::Int=30, N::Int=60, R::Float64=20.0)

Compute the elastic differential cross section at angle θ.
"""
function elastic_differential_cross_section(pot::OpticalPotential, E_cm::Float64, θ::Float64;
                                            l_max::Int=30, N::Int=60, R::Float64=20.0)
    μ = reduced_mass(pot)
    k = coulomb_k(μ, E_cm)
    η = coulomb_eta(pot.Z_proj, pot.Z_targ, μ, E_cm)

    sin_half = sin(θ/2)
    if abs(sin_half) < 1e-10
        return Inf
    end

    σ_0 = coulomb_sigma(0, η)
    f_c = -η / (2*k*sin_half^2) * exp(2im*σ_0 - 2im*η*log(sin_half))

    f_n = 0.0 + 0.0im
    cos_θ = cos(θ)

    for l in 0:l_max
        prob = ScatteringProblem(pot, E_cm, l; N=N, R=R)
        result = solve_scattering(prob)
        S_l = result.S_matrix

        σ_l = coulomb_sigma(l, η)
        P_l = legendre_P(l, cos_θ)

        f_n += (2*l + 1) * (S_l - 1) * exp(2im*σ_l) * P_l
    end
    f_n *= 1.0 / (2im*k)

    f_total = f_c + f_n
    return abs(f_total)^2 * 10.0
end

"""
    legendre_P(l::Int, x::Float64) -> Float64

Compute Legendre polynomial P_l(x) using recurrence relation.
"""
function legendre_P(l::Int, x::Float64)
    if l == 0
        return 1.0
    elseif l == 1
        return x
    else
        P_prev = 1.0
        P_curr = x
        for n in 2:l
            P_next = ((2*n - 1) * x * P_curr - (n - 1) * P_prev) / n
            P_prev = P_curr
            P_curr = P_next
        end
        return P_curr
    end
end
