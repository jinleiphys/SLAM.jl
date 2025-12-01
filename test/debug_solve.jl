# Detailed step-by-step debug of the scattering solver
using SLAM
using LinearAlgebra

function main()
    println("=" ^ 70)
    println("Step-by-step scattering solver debug")
    println("=" ^ 70)
    println()

    # Test parameters
    E_lab = 20.0
    A_proj = 1.0
    A_targ = 40.0
    E_cm = E_lab * A_targ / (A_proj + A_targ)

    pot = OpticalPotential(
        V_v=46.553, r_v=1.185, a_v=0.672,
        W_v=1.777, r_wv=1.185, a_wv=0.672,
        V_s=0.0, r_s=1.288, a_s=0.538,
        W_s=7.182, r_ws=1.288, a_ws=0.538,
        V_so=5.343, r_so=0.996, a_so=0.590,
        W_so=-0.110, r_wso=0.996, a_wso=0.590,
        r_c=0.0,
        Z_proj=0.0, A_proj=A_proj,
        Z_targ=0.0, A_targ=A_targ
    )

    l = 0
    j = 0.5
    N = 60  # Same as COLOSS
    R = 20.0

    # Constants
    HBARC = 197.327

    # Derived quantities
    μ = reduced_mass(pot)
    k = coulomb_k(μ, E_cm)
    η = coulomb_eta(pot.Z_proj, pot.Z_targ, μ, E_cm)
    coeff_kin = HBARC^2 / (2.0 * μ)

    println("Parameters:")
    println("  E_cm = $E_cm MeV")
    println("  μ = $μ MeV/c²")
    println("  k = $k fm⁻¹")
    println("  η = $η")
    println("  l = $l, j = $j")
    println("  N = $N, R = $R")
    println()

    # Initialize mesh
    mesh = init_legendre_mesh(N, R)

    # Compute Baye's T matrix
    T = baye_T_matrix(mesh)

    # Compute γ_s
    ρ_R = k * R
    γ_s = coulomb_gamma_s(l, η, ρ_R, k)

    println("Boundary condition parameter:")
    println("  ρ_R = kR = $ρ_R")
    println("  γ_s = $γ_s")
    println("  Expected for η=0: γ_s = ki = $(k*im)")
    println()

    # Build matrix and source
    M = zeros(ComplexF64, N, N)
    b = zeros(ComplexF64, N)

    println("Interior equations (rows 1 to $(N-1)):")
    println("-" ^ 70)
    for i in 1:N-1
        r_i = mesh.r[i]
        ρ_i = k * r_i

        V_i = evaluate_short_range(pot, r_i, l, j)
        U_i = V_i / coeff_kin

        F_l_i = coulomb_F(l, η, ρ_i)

        cent = Float64(l * (l + 1)) / r_i^2

        for jj in 1:N
            M[i, jj] = -T[i, jj]
            if i == jj
                M[i, jj] += k^2 - cent - U_i
            end
        end

        b[i] = U_i * F_l_i * sqrt(mesh.λ[i])

        println("Row $i: r = $(round(r_i, digits=4)), U = $(round(U_i, digits=4))")
        println("  k² - cent - U = $(round(k^2 - cent - U_i, digits=4))")
        println("  F_l(ρ) = $(round(F_l_i, digits=6)), b[$i] = $(round(b[i], digits=6))")
    end

    # Boundary condition row
    println()
    println("Boundary condition (row $N):")
    println("-" ^ 70)

    for jj in 1:N
        f_j_R = basis_function_at_R(mesh, jj)
        f_j_prime_R = basis_derivative_at_R(mesh, jj)
        M[N, jj] = f_j_prime_R - γ_s * f_j_R
        if jj <= 5
            println("  j=$jj: f_j(R) = $(round(f_j_R, digits=4)), f'_j(R) = $(round(f_j_prime_R, digits=4))")
            println("        M[$N,$jj] = f'_j(R) - γ_s*f_j(R) = $(round(M[N,jj], digits=4))")
        end
    end
    b[N] = 0.0

    println()
    println("Solving M * c = b...")
    println()

    # Solve
    c = M \ b

    println("Solution coefficients c:")
    for i in 1:N
        println("  c[$i] = $(round(c[i], digits=6))")
    end
    println()

    # Verify the solution
    residual = M * c - b
    println("Residual ||M*c - b|| = $(norm(residual))")
    println()

    # Extract φ(R)
    φ_R = sum(c[jj] * basis_function_at_R(mesh, jj) for jj in 1:N)

    # Also compute φ'(R)
    φ_prime_R = sum(c[jj] * basis_derivative_at_R(mesh, jj) for jj in 1:N)

    println("Wave function at boundary:")
    println("  φ(R) = $φ_R")
    println("  φ'(R) = $φ_prime_R")
    println("  φ'(R)/φ(R) = $(φ_prime_R / φ_R)")
    println("  Expected: φ'(R)/φ(R) = γ_s = $γ_s")
    println()

    # H+ at boundary
    H_plus_R = coulomb_H_plus(l, η, ρ_R)

    println("Coulomb functions at boundary:")
    println("  H⁺(kR) = $H_plus_R")
    println("  |H⁺(kR)| = $(abs(H_plus_R))")
    println()

    # S-matrix extraction
    f_l = φ_R / H_plus_R
    S_l = 1.0 + 2.0im * f_l

    println("S-matrix:")
    println("  f_l = φ(R) / H⁺(kR) = $f_l")
    println("  S_l = 1 + 2i*f_l = $S_l")
    println("  |S_l| = $(abs(S_l))")
    println()
    println("Expected COLOSS: S = (0.349464, 0.223937)")

    # Let me also try a different approach - using the asymptotic form
    println()
    println("=" ^ 70)
    println("Alternative analysis:")
    println("=" ^ 70)

    # For η=0, the outgoing wave is H⁺ = exp(ikr)
    # The scattered wave φ should have form: φ(r) → f_l * exp(ikr) as r → ∞
    # So f_l = φ(R) * exp(-ikR) for large R

    # Check: exp(-ikR)
    exp_minus_ikR = exp(-im * k * R)
    println("exp(-ikR) = $exp_minus_ikR")
    println("φ(R) * exp(-ikR) = $(φ_R * exp_minus_ikR)")
    println()

    # Another check: The S-matrix should satisfy |S| ≤ 1 for absorptive potential
    # and the imaginary part of the phase shift should be positive (absorption)

    δ_l = -0.5im * log(S_l)
    println("Phase shift: δ_l = $δ_l")
    println("  Re(δ_l) = $(real(δ_l))")
    println("  Im(δ_l) = $(imag(δ_l)) (should be > 0 for absorption)")
end

main()
