# Direct matrix element comparison for SLAM.jl debugging
# Same parameters, compare element by element

using SLAM
using Printf

function main()
    println("=" ^ 70)
    println("DIRECT MATRIX ELEMENT COMPARISON")
    println("=" ^ 70)

    # Test parameters
    # massp=1, masst=40, a1=0, a2=0 (so uses a13 = 40^(1/3))
    # zp=0, zt=0, elab=20

    E_lab = 20.0
    A_proj = 1.0
    A_targ = 40.0
    E_cm = E_lab * A_targ / (A_proj + A_targ)

    # Default convention: a1=0, a2=0 means use a13 = masst^(1/3) = 40^(1/3)
    a13_default = 40.0^(1/3)
    println("Default a13 = 40^(1/3) = $a13_default")

    # What SLAM uses with default (A1=0, A2=A_targ=40)
    a13_slam = 0.0^(1/3) + 40.0^(1/3)
    println("SLAM a13 = 0^(1/3) + 40^(1/3) = $a13_slam")
    println()

    # Create potential with default convention
    pot = OpticalPotential(
        V_v=46.553, r_v=1.185, a_v=0.672,
        W_v=1.777, r_wv=1.185, a_wv=0.672,
        V_s=0.0, r_s=1.288, a_s=0.538,
        W_s=7.182, r_ws=1.288, a_ws=0.538,
        V_so=5.343, r_so=0.996, a_so=0.590,
        W_so=-0.110, r_wso=0.996, a_wso=0.590,
        r_c=0.0,
        Z_proj=0.0, A_proj=A_proj,
        Z_targ=0.0, A_targ=A_targ,
        A1=0.0, A2=0.0  # Default convention
    )

    println("pot.A1 = $(pot.A1), pot.A2 = $(pot.A2)")
    println("Effective a13 = $(pot.A1^(1/3) + pot.A2^(1/3))")
    println()

    # Physical parameters
    HBARC = 197.327
    amu = 931.5
    m_proj = A_proj * amu
    m_targ = A_targ * amu
    μ = m_proj * m_targ / (m_proj + m_targ)
    k = sqrt(2 * μ * E_cm) / HBARC
    coeff_kin = HBARC^2 / (2.0 * μ)

    println("Physical parameters:")
    println("  E_cm = $E_cm MeV")
    println("  μ = $μ MeV/c²")
    println("  k = $k fm⁻¹")
    println("  coeff_kin = $coeff_kin MeV·fm²")
    println()

    # Same mesh
    N = 60
    R = 20.0
    l = 0
    j = 0.5

    mesh = init_legendre_mesh(N, R)
    T = baye_T_matrix(mesh)

    println("Mesh points (N=$N, R=$R):")
    for i in 1:N
        @printf("  r[%d] = %.10f fm\n", i, mesh.r[i])
    end
    println()

    println("T matrix (Baye's -d²/dr²):")
    for i in 1:N
        for jj in 1:N
            @printf(" %12.6f", T[i,jj])
        end
        println()
    end
    println()

    # Build matrix M and source b
    η = 0.0  # no Coulomb

    println("Matrix diagonal elements and source vector:")
    println("  i     r_i        U_i (real)      U_i (imag)      F_l         b_i (real)      b_i (imag)")
    println("-" ^ 100)

    for i in 1:N-1
        r_i = mesh.r[i]
        ρ_i = k * r_i

        V_i = evaluate_short_range(pot, r_i, l, j)
        U_i = V_i / coeff_kin

        F_l_i = coulomb_F(l, η, ρ_i)

        cent = Float64(l * (l + 1)) / r_i^2

        M_ii = -T[i,i] + k^2 - cent - U_i
        b_i = U_i * F_l_i * sqrt(mesh.λ[i])

        @printf("  %d  %10.6f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
                i, r_i, real(U_i), imag(U_i), F_l_i, real(b_i), imag(b_i))
    end
    println()

    # Now solve and get S-matrix
    prob = ScatteringProblem(pot, E_cm, l; j=j, N=N, R=R)
    result = solve_scattering(prob)

    println("S-matrix result:")
    println("  S = $(result.S_matrix)")
    println("  |S| = $(abs(result.S_matrix))")
end

main()
