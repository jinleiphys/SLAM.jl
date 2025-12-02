"""
    Numerov method for solving the radial Schrödinger equation.

Implements the Numerov algorithm for outward integration from r=0 to r=R,
with matching to asymptotic Coulomb functions to extract the S-matrix.

Reference: smoothie/scatt.f
"""

"""
    NumerovResult

Structure holding the results of a Numerov calculation.

# Fields
- `r::Vector{Float64}`: Radial grid points (fm)
- `ψ::Vector{ComplexF64}`: Radial wave function (unnormalized from Numerov)
- `ψ_normalized::Vector{ComplexF64}`: Normalized wave function ψ = (i/2)[H⁻ - S*H⁺]
- `S_matrix::ComplexF64`: S-matrix element
- `k::Float64`: Wave number (fm⁻¹)
- `η::Float64`: Sommerfeld parameter
- `σ_l::Float64`: Coulomb phase shift
"""
struct NumerovResult
    r::Vector{Float64}
    ψ::Vector{ComplexF64}
    ψ_normalized::Vector{ComplexF64}
    S_matrix::ComplexF64
    k::Float64
    η::Float64
    σ_l::Float64
end

"""
    solve_numerov(prob::ScatteringProblem; h::Float64=0.05) -> NumerovResult

Solve the scattering problem using the Numerov method.

The radial Schrödinger equation:
    [d²/dr² + k² - l(l+1)/r² - U(r)] ψ(r) = 0

where U(r) = 2μV(r)/ℏ².

Numerov algorithm (standard form):
    ψ_{n+1} = [2(1 - 5h²κ_n/12)ψ_n - (1 + h²κ_{n-1}/12)ψ_{n-1}] / (1 + h²κ_{n+1}/12)

where κ_n = k² - l(l+1)/r_n² - U(r_n)

# Arguments
- `prob::ScatteringProblem`: The scattering problem
- `h::Float64`: Step size (default 0.05 fm)

# Returns
- `NumerovResult`: Structure containing wave function and S-matrix
"""
function solve_numerov(prob::ScatteringProblem; h::Float64=0.05)
    pot = prob.pot
    E_cm = prob.E_cm
    l = prob.l
    j = prob.j
    R = prob.R

    # Compute physical parameters
    μ = reduced_mass(pot)
    k = coulomb_k(μ, E_cm)
    η = coulomb_eta(pot.Z_proj, pot.Z_targ, μ, E_cm)
    σ_l = coulomb_sigma(l, η)

    # Conversion factor
    coeff_kin = HBARC^2 / (2.0 * μ)

    # Check step size
    if k * h > 0.2
        @warn "Step size may be too large: k*h = $(k*h) > 0.2. Consider using h < $(0.2/k)"
    end

    # Setup grid: r[i] = i * h, i = 0, 1, ..., nr
    nr = Int(ceil(R / h))
    r = [(i * h) for i in 0:nr]

    # Starting point index: r0 = 2l * h (following smoothie/scatt.f convention)
    # r0_idx is the index where ψ = 0
    # For l=0: r0_idx = 0, so ψ[1] = 0 (at r=0)
    # For l>0: r0_idx = 2l, so ψ[r0_idx+1] = 0 (at r = 2l*h)
    r0_idx = 2 * l

    # Compute κ(r) = k² - l(l+1)/r² - U(r)
    # Use FULL potential (nuclear + Coulomb) for Numerov
    κ = zeros(ComplexF64, nr + 1)
    for i in 1:nr+1
        ri = r[i]
        if ri < 1e-10
            κ[i] = 0.0
        else
            V_full = evaluate_total_potential(pot, ri, l, j)
            U_full = V_full / coeff_kin
            κ[i] = k^2 - l * (l + 1) / ri^2 - U_full
        end
    end

    # Wave function array (1-indexed, r[1]=0, r[2]=h, ...)
    ψ = zeros(ComplexF64, nr + 1)

    # Boundary conditions (following smoothie/scatt.f):
    # ψ(r0) = 0  (at index r0_idx + 1)
    # ψ(r0 + h) = h^{l+1}  (at index r0_idx + 2)
    if r0_idx == 0
        # l = 0: ψ(0) = 0, ψ(h) = h
        ψ[1] = 0.0
        ψ[2] = h
    else
        ψ[r0_idx + 1] = 0.0
        ψ[r0_idx + 2] = h^(l + 1)
    end

    # For very high l, use smaller starting value
    if l > 100
        ψ[r0_idx + 2] = 1e-7
    end

    # Numerov algorithm coefficients
    h2 = h^2

    # For l=0, the first step needs special handling since κ[1] = κ(0) is undefined
    # Use simple Euler step for the first point if l=0
    start_idx = r0_idx + 2
    if r0_idx == 0
        # For l=0, do first step manually: ψ(2h) ≈ 2ψ(h) - h²κ(h)ψ(h)
        ψ[3] = 2 * ψ[2] - h2 * κ[2] * ψ[2]
        start_idx = 3
    end

    # Standard Numerov iteration (outward)
    # ψ_{n+1} = [2(1 - 5h²κ_n/12)ψ_n - (1 + h²κ_{n-1}/12)ψ_{n-1}] / (1 + h²κ_{n+1}/12)
    for i in start_idx:nr
        c_n = 1 - 5 * h2 * κ[i] / 12
        c_nm1 = 1 + h2 * κ[i - 1] / 12
        c_np1 = 1 + h2 * κ[i + 1] / 12
        ψ[i + 1] = (2 * c_n * ψ[i] - c_nm1 * ψ[i - 1]) / c_np1
    end

    # Matching to asymptotic form at r = R - 2h (to have points for derivative)
    # Use 5-point derivative formula: ψ'(R) ≈ (-ψ_{N} + 8ψ_{N-1} - 8ψ_{N-3} + ψ_{N-4}) / (12h)
    # Match at index nr-2
    i_match = nr - 2
    r_match = r[i_match + 1]  # r[nr-1] = (nr-2)*h

    # 5 points centered at i_match: indices i_match-2 to i_match+2 (1-indexed: i_match-1 to i_match+3)
    ψ_R = ψ[i_match + 1]
    ψp_R = (-ψ[i_match + 3] + 8*ψ[i_match + 2] - 8*ψ[i_match] + ψ[i_match - 1]) / (12 * h)

    # Get Coulomb functions at matching point
    ρ_match = k * r_match
    F_l, G_l, Fp_l, Gp_l = coulomb_FG(l, η, ρ_match)

    # H⁺ = G + iF, H⁻ = G - iF (outgoing and incoming Coulomb-Hankel functions)
    H_plus = G_l + im * F_l
    H_minus = G_l - im * F_l
    # Derivatives with respect to ρ
    Hp_plus = Gp_l + im * Fp_l
    Hp_minus = Gp_l - im * Fp_l

    # S-matrix from matching:
    # The normalized wave function asymptotically is:
    #   ψ_norm = (i/2) * [H⁻ - S * H⁺]
    #   ψ'_norm = (i/2) * k * [H'⁻ - S * H'⁺]
    # where derivatives are w.r.t. r, so dH/dr = k * dH/dρ
    #
    # Eliminating normalization:
    #   S = (ψ' * H⁻ - k * ψ * H'⁻) / (ψ' * H⁺ - k * ψ * H'⁺)
    S_matrix = (ψp_R * H_minus - k * ψ_R * Hp_minus) / (ψp_R * H_plus - k * ψ_R * Hp_plus)

    # Normalization factor N such that N * ψ = (i/2) * [H⁻ - S * H⁺]
    # At matching point: N * ψ_R = (i/2) * [H⁻ - S * H⁺]
    # So: N = (i/2) * [H⁻ - S * H⁺] / ψ_R
    N = (im / 2) * (H_minus - S_matrix * H_plus) / ψ_R

    # Normalized wave function: multiply by N and by e^{iσ_l}
    phase = exp(im * σ_l)
    ψ_normalized = phase * N * ψ

    return NumerovResult(r, ψ, ψ_normalized, S_matrix, k, η, σ_l)
end

"""
    solve_numerov_wavefunction(prob::ScatteringProblem; h::Float64=0.05) -> (NumerovResult, WaveFunctionResult)

Solve using both Numerov and Lagrange methods for comparison.
Returns the wave function in the same format for both methods.
"""
function compare_numerov_lagrange(prob::ScatteringProblem; h::Float64=0.05)
    # Solve with Numerov
    num_result = solve_numerov(prob; h=h)

    # Solve with Lagrange
    lag_result = solve_wavefunction(prob)

    return num_result, lag_result
end
