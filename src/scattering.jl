"""
    Main scattering solver using Lagrange-Legendre basis with exact D and T matrices.

Key insight: The D and T matrices act on expansion coefficients c_j,
NOT on function values φ(x_j)!

For x-regularized basis: f_j(r) = (r/r_j) * L_j(r) / √λ_j
Wave function: φ(r) = Σ_j c_j * f_j(r)
At mesh point r_j: f_j(r_j) = 1/√λ_j
So: c_j = φ(r_j) * √λ_j

Reference: D. Baye, Physics Reports 565 (2015) 1-107, Section 3.4.5
"""

# Physical constants (NIST 2014 values)
const HBARC = 197.3269718  # MeV·fm

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
    phase_shift::Float64          # Phase shift in degrees (real part)
    f_amplitude::ComplexF64
    cross_section::Float64
    converged::Bool
    l::Int                        # Angular momentum quantum number
    j::Float64                    # Total angular momentum
    k::Float64                    # Wave number (fm^-1)
    η::Float64                    # Sommerfeld parameter
end

"""
    solve_scattering(prob::ScatteringProblem) -> ScatteringResult

Solve the scattering problem using Lagrange mesh with coefficient-space formulation.

The radial Schrödinger equation:
[-ℏ²/(2μ) d²/dr² + ℏ²l(l+1)/(2μr²) + V(r)]φ = Eφ

Rewriting with U = 2μV/ℏ² and k² = 2μE/ℏ²:
[-d²/dr² + l(l+1)/r² + U(r) - k²]φ = 0

For Coulomb scattering, we decompose the wave function:
ψ(r) = F_l(η, kr) + ψ_sc(r)

The scattered wave satisfies:
[-d²/dr² + l(l+1)/r² + U_full - k²]ψ_sc = U_short * F_l(η, kr)

Key points:
- Matrix M uses FULL potential: U_full = U_nuc + U_coul_finite
- Source term b uses SHORT-RANGE potential: U_short = U_nuc + U_coul_finite - U_coul_point
- Phase factors e^{iσ_l} cancel out in this formulation
- S-matrix: S = 1 + 2ik*f_l (no extra phase factor)

Reference: D. Baye, Physics Reports 565 (2015) 1-107
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

        # Matrix uses FULL potential: V_full = V_nuc + V_coul_finite
        # Source uses SHORT-RANGE potential: V_short = V_nuc + V_coul_finite - V_coul_point
        V_full_i = evaluate_total_potential(pot, r_i, l, j)
        V_short_i = evaluate_short_range(pot, r_i, l, j)
        U_full_i = V_full_i / coeff_kin
        U_short_i = V_short_i / coeff_kin

        # Evaluate Coulomb function F_l at ρ_i
        F_l_i = coulomb_F(l, η, ρ_i)

        # Centrifugal term
        cent = Float64(l * (l + 1)) / r_i^2

        # Build matrix row using FULL potential
        # M = -T + (k² - l(l+1)/r² - U_full) on diagonal
        for jj in 1:N
            M[i, jj] = -T[i, jj]
            if i == jj
                M[i, jj] += k^2 - cent - U_full_i
            end
        end

        # Source term uses SHORT-RANGE potential: U_short * F_l * √λ
        # Note: No exp(iσ_l) factor - it cancels in this formulation
        b[i] = U_short_i * F_l_i * sqrt(mesh.λ[i])
    end

    # Boundary condition at r = R (row N)
    # φ'(R) - γ_s * φ(R) = 0
    # Using analytical formulas:
    #   f̂_j(1) = (-1)^{N-j} / √[x_j(1-x_j)]
    #   df̂_j/dx|_{x=1} = (-1)^{N-j}/√[x_j(1-x_j)] * [N(N+1) - x_j/(1-x_j)]
    # Converting to r-coordinates: f_j(R) = f̂_j(1), f_j'(R) = (1/R) df̂_j/dx
    for jj in 1:N
        f_j_R = fhat_at_boundary(mesh, jj)
        f_j_prime_R = dfhat_dx_at_boundary(mesh, jj) / R
        M[N, jj] = f_j_prime_R - γ_s * f_j_R
    end
    b[N] = 0.0

    # Solve the linear system
    c = M \ b

    # Extract scattering amplitude using direct boundary matching
    # φ_sc(R) = Σ_j c_j * f_j(R)
    φ_R = sum(c[jj] * basis_function_at_R(mesh, jj) for jj in 1:N)

    # H⁺_l(kR) = G_l + i*F_l (outgoing Coulomb-Hankel function)
    H_plus_R = coulomb_H_plus(l, η, ρ_R)

    # Scattering amplitude: f_l = φ_sc(R) / [k * H⁺_l(kR)]
    # Note: No exp(iσ_l) factor needed - phases cancel in this formulation
    f_l = φ_R / (k * H_plus_R)

    # S-matrix: S_l = 1 + 2ik*f_l
    # The phase factors are absorbed in the definition of f_l
    S_l = 1.0 + 2.0im * k * f_l

    # Phase shift: S_l = exp(2iδ_l)
    # δ_l = arg(S_l)/2 in radians, convert to degrees
    δ_l_rad = angle(S_l) / 2
    δ_l_deg = rad2deg(δ_l_rad)

    # Partial wave cross section
    # σ_reaction = (π/k²) * (2J+1)/(2S+1) * (1 - |S_l|²) * 10 (fm² to mb)
    s_spin = 0.5  # spin-1/2 projectile
    σ_l = π / k^2 / (2*s_spin + 1) * (2*j + 1) * (1.0 - abs(S_l)^2) * 10.0

    return ScatteringResult(S_l, δ_l_deg, f_l, σ_l, true, l, j, k, η)
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
                            l_max::Int=20, N::Int=60, R::Float64=20.0,
                            S_threshold::Float64=0.999)

Solve scattering for all partial waves until |S_l| > S_threshold or l > l_max.

Returns a Vector of ScatteringResult for all computed partial waves.
"""
function solve_all_partial_waves(pot::OpticalPotential, E_cm::Float64;
                                  l_max::Int=20, N::Int=60, R::Float64=20.0,
                                  S_threshold::Float64=0.999)
    results = ScatteringResult[]

    for l in 0:l_max
        if l == 0
            j_values = [0.5]
        else
            j_values = [Float64(l) - 0.5, Float64(l) + 0.5]
        end

        all_converged = true
        for j in j_values
            prob = ScatteringProblem(pot, E_cm, l; j=j, N=N, R=R)
            result = solve_scattering(prob)
            push!(results, result)

            # Check if S-matrix magnitude is close to 1 (no absorption)
            if abs(result.S_matrix) < S_threshold
                all_converged = false
            end
        end

        # Early termination if all j values at this l have |S| ≈ 1
        if all_converged && l > 5
            break
        end
    end

    return results
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

"""
    WaveFunctionResult

Structure holding the wave function calculation results.

# Fields
- `r::Vector{Float64}`: Radial grid points (fm)
- `ψ_sc::Vector{ComplexF64}`: Scattered wave function ψ^{sc}(r)
- `ψ_total::Vector{ComplexF64}`: Total wave function ψ(r) = e^{iσ_l}[F_l + ψ^{sc}]
- `F_l::Vector{Float64}`: Regular Coulomb function F_l(η, kr)
- `c::Vector{ComplexF64}`: Expansion coefficients in Lagrange basis
- `σ_l::Float64`: Coulomb phase shift σ_l
- `S_matrix::ComplexF64`: S-matrix element
"""
struct WaveFunctionResult
    r::Vector{Float64}
    ψ_sc::Vector{ComplexF64}
    ψ_total::Vector{ComplexF64}
    F_l::Vector{Float64}
    c::Vector{ComplexF64}
    σ_l::Float64
    S_matrix::ComplexF64
end

"""
    solve_wavefunction(prob::ScatteringProblem) -> WaveFunctionResult

Solve the scattering problem and return the full wave function.

The total wave function is:
    ψ(r) = e^{iσ_l} [F_l(η, kr) + ψ^{sc}(r)]

where:
- σ_l is the Coulomb phase shift
- F_l(η, kr) is the regular Coulomb function (incident wave)
- ψ^{sc}(r) is the scattered wave

Note: In our formulation, ψ^{sc} is computed WITHOUT the e^{iσ_l} factor,
so the total wave function requires multiplying by e^{iσ_l}.

# Arguments
- `prob::ScatteringProblem`: The scattering problem specification

# Returns
- `WaveFunctionResult`: Structure containing wave functions and coefficients
"""
function solve_wavefunction(prob::ScatteringProblem)
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

    # Coulomb phase shift
    σ_l = coulomb_sigma(l, η)

    # Conversion factor: U = 2μV/ℏ² (converts V in MeV to U in fm⁻²)
    coeff_kin = HBARC^2 / (2.0 * μ)

    # Initialize mesh
    mesh = init_legendre_mesh(N, R)

    # Compute Baye's T matrix
    T = baye_T_matrix(mesh)

    # Compute boundary condition parameter γ_s
    ρ_R = k * R
    γ_s = coulomb_gamma_s(l, η, ρ_R, k)

    # Build the matrix equation M * c = b
    M = zeros(ComplexF64, N, N)
    b = zeros(ComplexF64, N)

    # Store F_l values at mesh points
    F_l_mesh = zeros(Float64, N)

    # Interior points (i = 1 to N-1)
    for i in 1:N-1
        r_i = mesh.r[i]
        ρ_i = k * r_i

        # Potentials
        V_full_i = evaluate_total_potential(pot, r_i, l, j)
        V_short_i = evaluate_short_range(pot, r_i, l, j)
        U_full_i = V_full_i / coeff_kin
        U_short_i = V_short_i / coeff_kin

        # Coulomb function
        F_l_mesh[i] = coulomb_F(l, η, ρ_i)

        # Centrifugal term
        cent = Float64(l * (l + 1)) / r_i^2

        # Build matrix row
        for jj in 1:N
            M[i, jj] = -T[i, jj]
            if i == jj
                M[i, jj] += k^2 - cent - U_full_i
            end
        end

        # Source term
        b[i] = U_short_i * F_l_mesh[i] * sqrt(mesh.λ[i])
    end

    # Last mesh point F_l
    F_l_mesh[N] = coulomb_F(l, η, k * mesh.r[N])

    # Boundary condition at r = R (using analytical formulas)
    for jj in 1:N
        f_j_R = fhat_at_boundary(mesh, jj)
        f_j_prime_R = dfhat_dx_at_boundary(mesh, jj) / R
        M[N, jj] = f_j_prime_R - γ_s * f_j_R
    end
    b[N] = 0.0

    # Solve the linear system
    c = M \ b

    # Compute wave functions at mesh points
    ψ_sc = zeros(ComplexF64, N)
    ψ_total = zeros(ComplexF64, N)

    # Phase factor e^{iσ_l}
    phase = exp(im * σ_l)

    for i in 1:N
        # ψ^{sc}(r_i) = c_i / √λ_i (from c_i = ψ^{sc}(r_i) * √λ_i)
        ψ_sc[i] = c[i] / sqrt(mesh.λ[i])

        # Total wave function: ψ(r) = e^{iσ_l} [F_l + ψ^{sc}]
        ψ_total[i] = phase * (F_l_mesh[i] + ψ_sc[i])
    end

    # Compute S-matrix (using analytical formula for f̂_j(1))
    φ_R = sum(c[jj] * fhat_at_boundary(mesh, jj) for jj in 1:N)
    H_plus_R = coulomb_H_plus(l, η, ρ_R)
    f_l = φ_R / (k * H_plus_R)
    S_matrix = 1.0 + 2.0im * k * f_l

    return WaveFunctionResult(copy(mesh.r), ψ_sc, ψ_total, F_l_mesh, c, σ_l, S_matrix)
end

"""
    evaluate_wavefunction(wf::WaveFunctionResult, r::Float64) -> ComplexF64

Evaluate the total wave function at a single radial point using the stored coefficients.
Uses Lagrange interpolation through the mesh points.

Returns ψ_total(r) = e^{iσ_l} [F_l(η, kr) + ψ^{sc}(r)]
"""
function evaluate_wavefunction(wf::WaveFunctionResult, r::Float64)
    # Simple linear interpolation for points within the mesh
    N = length(wf.r)

    if r <= 0
        return 0.0 + 0.0im
    end

    R = wf.r[end] + (wf.r[end] - wf.r[end-1])  # Approximate R

    if r > R
        return wf.ψ_total[end]  # Return last value
    end

    # Find bracketing indices
    idx = searchsortedfirst(wf.r, r)
    if idx == 1
        return wf.ψ_total[1]
    elseif idx > N
        return wf.ψ_total[N]
    end

    # Linear interpolation
    r1, r2 = wf.r[idx-1], wf.r[idx]
    ψ1, ψ2 = wf.ψ_total[idx-1], wf.ψ_total[idx]
    t = (r - r1) / (r2 - r1)

    return ψ1 + t * (ψ2 - ψ1)
end

"""
    evaluate_wavefunction(wf::WaveFunctionResult, r_eval::Vector{Float64},
                          prob::ScatteringProblem) -> (Vector{ComplexF64}, Vector{ComplexF64})

Evaluate the wave function at arbitrary radial points using the Lagrange basis.

Returns (ψ_sc, ψ_total) at the evaluation points.
"""
function evaluate_wavefunction(wf::WaveFunctionResult, r_eval::Vector{Float64},
                               prob::ScatteringProblem)
    pot = prob.pot
    l = prob.l
    N = prob.N
    R = prob.R

    μ = reduced_mass(pot)
    k = coulomb_k(μ, prob.E_cm)
    η = coulomb_eta(pot.Z_proj, pot.Z_targ, μ, prob.E_cm)

    mesh = init_legendre_mesh(N, R)

    n_eval = length(r_eval)
    ψ_sc_eval = zeros(ComplexF64, n_eval)
    ψ_total_eval = zeros(ComplexF64, n_eval)

    phase = exp(im * wf.σ_l)

    for (idx, r) in enumerate(r_eval)
        if r <= 0 || r > R
            continue
        end

        # ψ^{sc}(r) = Σ_j c_j * f_j(r)
        # where f_j(r) = (r/r_j) * L_j(r) / √λ_j
        for j in 1:N
            L_j = lagrange_polynomial(mesh, j, r)
            f_j = (r / mesh.r[j]) * L_j / sqrt(mesh.λ[j])
            ψ_sc_eval[idx] += wf.c[j] * f_j
        end

        # F_l at this point
        F_l_r = coulomb_F(l, η, k * r)

        # Total: ψ = e^{iσ_l} [F_l + ψ^{sc}]
        ψ_total_eval[idx] = phase * (F_l_r + ψ_sc_eval[idx])
    end

    return ψ_sc_eval, ψ_total_eval
end
