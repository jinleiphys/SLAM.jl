# Debug script to compare intermediate values in SLAM.jl
using SLAM
using Printf

function main()
    println("=" ^ 70)
    println("Debug: Intermediate Value Comparison")
    println("=" ^ 70)

    # Test parameters (no Coulomb case)
    E_lab = 20.0
    A_proj = 1.0
    A_targ = 40.0
    E_cm = E_lab * A_targ / (A_proj + A_targ)

    # Create optical potential
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

    # Constants
    HBARC = 197.327

    # Compute derived quantities
    μ = reduced_mass(pot)
    k = coulomb_k(μ, E_cm)
    η = coulomb_eta(pot.Z_proj, pot.Z_targ, μ, E_cm)
    coeff_kin = HBARC^2 / (2.0 * μ)

    println("Physical parameters:")
    @printf("  E_cm = %.6f MeV\n", E_cm)
    @printf("  μ = %.6f MeV/c²\n", μ)
    @printf("  k = %.6f fm⁻¹\n", k)
    @printf("  η = %.6f\n", η)
    @printf("  coeff_kin = %.6f MeV·fm²\n", coeff_kin)
    println()

    # Initialize mesh
    N = 60
    R = 20.0
    mesh = init_legendre_mesh(N, R)

    println("First few mesh points:")
    for j in 1:5
        @printf("  r[%d] = %.10f, x[%d] = %.10f, λ[%d] = %.10f\n",
                j, mesh.r[j], j, mesh.x[j], j, mesh.λ[j])
    end
    println()

    # Test Baye T matrix
    T = baye_T_matrix(mesh)
    println("Baye T matrix (first 3x3 block):")
    for i in 1:3
        @printf("  ")
        for j in 1:3
            @printf("%12.6f ", T[i, j])
        end
        println()
    end
    println()

    # Test potential at a few points
    l = 0
    j_qn = 0.5

    println("Potential values at first few mesh points (l=$l, j=$j_qn):")
    println("  i       r          V_nuc (MeV)                      U (fm⁻²)")
    for i in 1:5
        r_i = mesh.r[i]
        V_nuc = evaluate_potential(pot, r_i, l, j_qn)
        U_i = V_nuc / coeff_kin
        @printf("  %d  %.6f  (%10.4f, %10.4f)  (%12.6f, %12.6f)\n",
                i, r_i, real(V_nuc), imag(V_nuc), real(U_i), imag(U_i))
    end
    println()

    # Test Coulomb functions
    println("Coulomb functions at first few points:")
    for i in 1:5
        ρ_i = k * mesh.r[i]
        F_l = coulomb_F(l, η, ρ_i)
        @printf("  ρ[%d] = %.6f, F_%d = %.10f\n", i, ρ_i, l, F_l)
    end
    println()

    # Test boundary condition parameter
    ρ_R = k * R
    γ_s = coulomb_gamma_s(l, η, ρ_R)
    @printf("Boundary condition: γ_s = (%.6f, %.6f)\n", real(γ_s), imag(γ_s))
    println()

    # Test basis functions at R
    println("Basis functions at R:")
    for jj in 1:5
        f_j_R = basis_function_at_R(mesh, jj)
        f_j_prime_R = basis_derivative_at_R(mesh, jj)
        @printf("  f_%d(R) = %.10f, f'_%d(R) = %.10f\n", jj, f_j_R, jj, f_j_prime_R)
    end
    println()

    # Build matrix and source vector (first few rows)
    M = zeros(ComplexF64, N, N)
    b = zeros(ComplexF64, N)

    for i in 1:N-1
        r_i = mesh.r[i]
        ρ_i = k * r_i

        V_i = evaluate_short_range(pot, r_i, l, j_qn)
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
    end

    # Boundary condition
    for jj in 1:N
        f_j_R = basis_function_at_R(mesh, jj)
        f_j_prime_R = basis_derivative_at_R(mesh, jj)
        M[N, jj] = f_j_prime_R - γ_s * f_j_R
    end
    b[N] = 0.0

    println("Matrix M diagonal elements (first 5):")
    for i in 1:5
        @printf("  M[%d,%d] = (%12.6f, %12.6f)\n", i, i, real(M[i, i]), imag(M[i, i]))
    end
    println()

    println("Source vector b (first 5):")
    for i in 1:5
        @printf("  b[%d] = (%12.6f, %12.6f)\n", i, real(b[i]), imag(b[i]))
    end
    println()

    # Solve
    c = M \ b

    println("Solution coefficients c (first 5):")
    for i in 1:5
        @printf("  c[%d] = (%12.6f, %12.6f)\n", i, real(c[i]), imag(c[i]))
    end
    println()

    # Extract phi(R) and S-matrix
    φ_R = sum(c[jj] * basis_function_at_R(mesh, jj) for jj in 1:N)
    H_plus_R = coulomb_H_plus(l, η, ρ_R)

    f_l = φ_R / H_plus_R
    S_l = 1.0 + 2.0im * f_l

    @printf("φ(R) = (%.10f, %.10f)\n", real(φ_R), imag(φ_R))
    @printf("H⁺(kR) = (%.10f, %.10f)\n", real(H_plus_R), imag(H_plus_R))
    @printf("f_l = (%.10f, %.10f)\n", real(f_l), imag(f_l))
    @printf("S_l = (%.10f, %.10f), |S| = %.10f\n", real(S_l), imag(S_l), abs(S_l))
end

main()
