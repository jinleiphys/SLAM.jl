"""
Compare Numerov and SLAM methods for p + 12C scattering using KD02 potential.

This script:
1. Sets up p + 12C scattering with KD02 global optical potential
2. Compares wave functions from Numerov and SLAM methods
3. Compares S-matrix elements for multiple partial waves
4. Outputs data to files for Python plotting
"""

using Printf
using LinearAlgebra
using DelimitedFiles

# Add SLAM.jl to path and load
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using SLAM

# Load KD02 potential module
include(joinpath(@__DIR__, "..", "src", "kd02.jl"))
using .KD02

# ============================================================================
# Setup: p + 12C scattering parameters
# ============================================================================

const Z_proj = 1.0    # proton
const A_proj = 1.0
const Z_targ = 6.0    # 12C
const A_targ = 12.0
const Elab = 30.0     # MeV (laboratory energy)

# Convert Elab to Ecm
function lab_to_cm(Elab, A_proj, A_targ)
    return Elab * A_targ / (A_proj + A_targ)
end

const E_cm = lab_to_cm(Elab, A_proj, A_targ)

println("="^60)
println("p + 12C Scattering: Numerov vs SLAM Comparison")
println("="^60)
@printf("Laboratory Energy: %.2f MeV\n", Elab)
@printf("Center-of-mass Energy: %.4f MeV\n", E_cm)
println()

# ============================================================================
# Get KD02 potential parameters
# ============================================================================

kd02 = kd02_potential(2, Int(Z_targ), Int(A_targ), Elab, verbose=true)

# Create OpticalPotential from KD02 parameters (spin=0, so no spin-orbit)
pot = OpticalPotential(
    V_v = kd02.v,      r_v = kd02.rv,   a_v = kd02.av,
    W_v = kd02.w,      r_wv = kd02.rw,  a_wv = kd02.aw,
    V_s = kd02.vd,     r_s = kd02.rvd,  a_s = kd02.avd,
    W_s = kd02.wd,     r_ws = kd02.rwd, a_ws = kd02.awd,
    V_so = 0.0,        r_so = kd02.rvso, a_so = kd02.avso,  # spin=0
    W_so = 0.0,        r_wso = kd02.rwso, a_wso = kd02.awso,
    r_c = kd02.rc,
    Z_proj = Z_proj, A_proj = A_proj,
    Z_targ = Z_targ, A_targ = A_targ
)

println("\nOptical Potential Summary:")
@printf("  Real volume:    V = %.3f MeV, r = %.3f fm, a = %.3f fm\n", pot.V_v, pot.r_v, pot.a_v)
@printf("  Imag volume:    W = %.3f MeV, r = %.3f fm, a = %.3f fm\n", pot.W_v, pot.r_wv, pot.a_wv)
@printf("  Imag surface:   Wd = %.3f MeV, r = %.3f fm, a = %.3f fm\n", pot.W_s, pot.r_ws, pot.a_ws)
@printf("  Coulomb:        rc = %.3f fm\n", pot.r_c)
println()

# ============================================================================
# Scattering calculations
# ============================================================================

# Parameters
const N_mesh = 80      # Number of Lagrange mesh points
const R_max = 25.0     # Maximum radius (fm)
const h_numerov = 0.02 # Numerov step size (fm)
const l_max = 10       # Maximum partial wave

# Compute reduced mass and wave number
mu = reduced_mass(pot)
k = coulomb_k(mu, E_cm)
eta = coulomb_eta(Z_proj, Z_targ, mu, E_cm)

@printf("Physical parameters:\n")
@printf("  Reduced mass: mu = %.4f MeV/c^2\n", mu)
@printf("  Wave number:  k = %.6f fm^-1\n", k)
@printf("  Sommerfeld:   eta = %.6f\n", eta)
println()

# ============================================================================
# Compare S-matrix for multiple partial waves
# ============================================================================

println("="^60)
println("S-matrix Comparison (Numerov vs SLAM)")
println("="^60)
println()
@printf("%3s  %18s  %18s  %18s  %10s\n",
        "l", "S_Numerov", "S_SLAM", "Delta|S|", "Delta_arg(S)")
println("-"^80)

S_numerov = ComplexF64[]
S_slam = ComplexF64[]
l_values = Int[]

for l in 0:l_max
    j = Float64(l)
    prob = ScatteringProblem(pot, E_cm, l; j=j, N=N_mesh, R=R_max)
    num_result = solve_numerov(prob; h=h_numerov)
    slam_result = solve_scattering(prob)

    push!(S_numerov, num_result.S_matrix)
    push!(S_slam, slam_result.S_matrix)
    push!(l_values, l)

    delta_mod = abs(abs(num_result.S_matrix) - abs(slam_result.S_matrix))
    delta_arg = abs(angle(num_result.S_matrix) - angle(slam_result.S_matrix)) * 180 / pi

    @printf("%3d  %8.5f%+8.5fi  %8.5f%+8.5fi  %10.2e  %8.4f deg\n",
            l,
            real(num_result.S_matrix), imag(num_result.S_matrix),
            real(slam_result.S_matrix), imag(slam_result.S_matrix),
            delta_mod, delta_arg)
end
println()

# ============================================================================
# Compare wave functions for selected partial waves
# ============================================================================

println("="^60)
println("Wave Function Comparison")
println("="^60)

l_compare = [0, 2, 5]
wf_data = Dict()

for l in l_compare
    j = Float64(l)
    prob = ScatteringProblem(pot, E_cm, l; j=j, N=N_mesh, R=R_max)
    num_result = solve_numerov(prob; h=h_numerov)
    wf_result = solve_wavefunction(prob)

    wf_data[l] = (numerov = num_result, slam = wf_result, prob = prob)

    println("\nPartial wave l = $l:")
    @printf("%8s  %24s  %24s  %12s\n", "r (fm)", "psi_Numerov", "psi_SLAM", "Delta(rel)")
    println("-"^72)

    r_compare = [2.0, 5.0, 10.0, 15.0, 20.0]
    for r in r_compare
        if r > R_max continue end
        idx_num = argmin(abs.(num_result.r .- r))
        psi_num = num_result.ψ_normalized[idx_num]
        r_num = num_result.r[idx_num]
        idx_slam = argmin(abs.(wf_result.r .- r))
        psi_slam = wf_result.ψ_total[idx_slam]
        delta_rel = abs(psi_num - psi_slam) / max(abs(psi_num), abs(psi_slam), 1e-10)
        @printf("%8.2f  %11.5f%+11.5fi  %11.5f%+11.5fi  %12.2e\n",
                r_num, real(psi_num), imag(psi_num), real(psi_slam), imag(psi_slam), delta_rel)
    end
end

# ============================================================================
# Output data to files for Python plotting
# ============================================================================

println("\n" * "="^60)
println("Outputting data for Python plotting")
println("="^60)

data_dir = joinpath(@__DIR__, "..", "figures", "data")
mkpath(data_dir)

# S-matrix data
open(joinpath(data_dir, "smatrix.dat"), "w") do f
    println(f, "# p + 12C scattering, E_lab = $(Elab) MeV")
    println(f, "# l  Re(S_num)  Im(S_num)  Re(S_slam)  Im(S_slam)")
    for i in eachindex(l_values)
        @printf(f, "%d  %.10f  %.10f  %.10f  %.10f\n",
                l_values[i], real(S_numerov[i]), imag(S_numerov[i]),
                real(S_slam[i]), imag(S_slam[i]))
    end
end
println("Saved: figures/data/smatrix.dat")

# Wave function data for each l
r_max_plot = 15.0
for l in l_compare
    data = wf_data[l]
    num = data.numerov
    slam = data.slam

    # Numerov data
    idx_num = num.r .<= r_max_plot
    open(joinpath(data_dir, "wf_numerov_l$(l).dat"), "w") do f
        println(f, "# Numerov wave function for l=$l")
        println(f, "# r  Re(psi)  Im(psi)")
        for i in findall(idx_num)
            @printf(f, "%.6f  %.10f  %.10f\n",
                    num.r[i], real(num.ψ_normalized[i]), imag(num.ψ_normalized[i]))
        end
    end

    # SLAM data
    idx_slam = slam.r .<= r_max_plot
    open(joinpath(data_dir, "wf_slam_l$(l).dat"), "w") do f
        println(f, "# SLAM wave function for l=$l")
        println(f, "# r  Re(psi)  Im(psi)")
        for i in findall(idx_slam)
            @printf(f, "%.6f  %.10f  %.10f\n",
                    slam.r[i], real(slam.ψ_total[i]), imag(slam.ψ_total[i]))
        end
    end
    println("Saved: figures/data/wf_numerov_l$(l).dat and wf_slam_l$(l).dat")
end

println("\nAll data files saved. Run plot_p12C.py to generate figures.")
println("="^60)
