#!/usr/bin/env julia
"""
Compare source term integration methods:
1. LAGRANGE_GAUSS: Original Lagrange-Gauss quadrature at mesh points
2. FINE_GRID: Finer grid with Lagrange basis interpolation

This script compares the S-matrix results from both methods to assess
the accuracy improvement from using finer grid integration.
"""

using SLAM
using Printf

println("="^70)
println("Comparison of Source Term Integration Methods")
println("="^70)

# Define potential (p + 12C at 30.3 MeV)
pot = OpticalPotential(
    V_v=53.30, r_v=1.25, a_v=0.65,
    W_v=0.0, r_wv=1.25, a_wv=0.47,
    W_s=9.31, r_ws=1.25, a_ws=0.47,  # Surface absorption
    V_so=6.20, r_so=1.25, a_so=0.65,
    W_so=0.0, r_wso=1.25, a_wso=0.65,
    r_c=1.25,
    Z_proj=1.0, A_proj=1.0,
    Z_targ=6.0, A_targ=12.0
)

E_cm = 30.3
R = 15.0

println("\nPotential: p + 12C at E_cm = $E_cm MeV")
println("Matching radius R = $R fm")

# Test different N values (Lagrange mesh points)
N_values = [20, 30, 40, 50, 60]
n_fine_values = [100, 200, 300, 500]

println("\n" * "="^70)
println("Part 1: Compare methods for different N (Lagrange mesh points)")
println("="^70)

println("\n" * "-"^70)
@printf("%-5s | %-25s | %-25s | %-12s\n", "N", "S (Lagrange-Gauss)", "S (Fine Grid, n=200)", "Diff |S|")
println("-"^70)

for N in N_values
    prob = ScatteringProblem(pot, E_cm, 0; N=N, R=R)

    result_lg = solve_scattering(prob; source_method=LAGRANGE_GAUSS)
    result_fg = solve_scattering(prob; source_method=FINE_GRID, n_fine=200)

    S_lg = result_lg.S_matrix
    S_fg = result_fg.S_matrix
    diff = abs(abs(S_lg) - abs(S_fg))

    @printf("%-5d | %10.6f + %10.6fi | %10.6f + %10.6fi | %12.2e\n",
            N, real(S_lg), imag(S_lg), real(S_fg), imag(S_fg), diff)
end

println("\n" * "="^70)
println("Part 2: Effect of n_fine (fine grid points) for N=40")
println("="^70)

N = 40
prob = ScatteringProblem(pot, E_cm, 0; N=N, R=R)
result_lg = solve_scattering(prob; source_method=LAGRANGE_GAUSS)
S_ref = result_lg.S_matrix

println("\nReference (Lagrange-Gauss): S = $(real(S_ref)) + $(imag(S_ref))i")
println("\n" * "-"^70)
@printf("%-10s | %-25s | %-12s | %-12s\n", "n_fine", "S (Fine Grid)", "|S| diff", "Phase diff")
println("-"^70)

for n_fine in n_fine_values
    result_fg = solve_scattering(prob; source_method=FINE_GRID, n_fine=n_fine)
    S_fg = result_fg.S_matrix

    diff_abs = abs(abs(S_fg) - abs(S_ref))
    diff_phase = abs(angle(S_fg) - angle(S_ref)) * 180 / π

    @printf("%-10d | %10.6f + %10.6fi | %12.2e | %10.4f°\n",
            n_fine, real(S_fg), imag(S_fg), diff_abs, diff_phase)
end

println("\n" * "="^70)
println("Part 3: Comparison across all partial waves (l=0 to 10)")
println("="^70)

N = 40
println("\nN = $N, n_fine = 200")
println("\n" * "-"^90)
@printf("%-3s | %-25s | %-25s | %-12s | %-12s\n",
        "l", "S (Lagrange-Gauss)", "S (Fine Grid)", "|S| diff", "δ diff (°)")
println("-"^90)

for l in 0:10
    local prob = ScatteringProblem(pot, E_cm, l; N=N, R=R)

    local result_lg = solve_scattering(prob; source_method=LAGRANGE_GAUSS)
    local result_fg = solve_scattering(prob; source_method=FINE_GRID, n_fine=200)

    S_lg = result_lg.S_matrix
    S_fg = result_fg.S_matrix

    diff_abs = abs(abs(S_lg) - abs(S_fg))
    diff_phase = abs(result_lg.phase_shift - result_fg.phase_shift)

    @printf("%-3d | %10.6f + %10.6fi | %10.6f + %10.6fi | %12.2e | %10.4f\n",
            l, real(S_lg), imag(S_lg), real(S_fg), imag(S_fg), diff_abs, diff_phase)
end

println("\n" * "="^70)
println("Part 4: Source term vector comparison (l=0, N=40)")
println("="^70)

N = 40
prob = ScatteringProblem(pot, E_cm, 0; N=N, R=R)

# Get mesh and parameters for source term calculation
μ = SLAM.reduced_mass(pot)
k = SLAM.coulomb_k(μ, E_cm)
η = SLAM.coulomb_eta(pot.Z_proj, pot.Z_targ, μ, E_cm)
coeff_kin = SLAM.HBARC^2 / (2.0 * μ)
mesh = SLAM.init_legendre_mesh(N, R)

# Compute source terms with both methods
b_lg = zeros(ComplexF64, N-1)
for i in 1:N-1
    r_i = mesh.r[i]
    ρ_i = k * r_i
    V_short_i = SLAM.evaluate_short_range(pot, r_i, 0, 0.5)
    U_short_i = V_short_i / coeff_kin
    F_l_i = SLAM.coulomb_F(0, η, ρ_i)
    b_lg[i] = U_short_i * F_l_i * sqrt(R * mesh.λ[i])
end

b_fg = compute_source_term_fine(mesh, pot, 0, 0.5, k, η, coeff_kin; n_fine=200)

println("\nSource term b[j] comparison (first 10 elements):")
println("-"^70)
@printf("%-5s | %-12s | %-20s | %-20s | %-12s\n",
        "j", "r[j] (fm)", "b_LG", "b_FG", "Rel. diff")
println("-"^70)

for j in 1:min(10, N-1)
    rel_diff = abs(b_lg[j] - b_fg[j]) / max(abs(b_lg[j]), 1e-15)
    @printf("%-5d | %12.4f | %20.10f | %20.10f | %12.2e\n",
            j, mesh.r[j], real(b_lg[j]), real(b_fg[j]), rel_diff)
end

println("\n" * "="^70)
println("Summary")
println("="^70)
println("""
The comparison shows:
1. Both methods give similar results for large N (well-resolved potential)
2. Fine grid integration may improve accuracy when potential varies rapidly
3. Relative differences are typically small (< 10^-4) for reasonable N values
4. The Lagrange-Gauss method is faster but may miss rapid potential variations
""")
