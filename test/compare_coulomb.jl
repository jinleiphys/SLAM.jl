# Test Coulomb treatment in SLAM.jl
# Compare the short-range potential calculation with expected behavior

using SLAM
using Printf

function main()
    println("=" ^ 70)
    println("COULOMB TREATMENT TEST")
    println("=" ^ 70)
    println()

    # Test parameters: p + 40Ca at 30 MeV
    A_proj = 1.0
    A_targ = 40.0
    Z_proj = 1.0
    Z_targ = 20.0
    E_lab = 30.0
    E_cm = E_lab * A_targ / (A_proj + A_targ)

    println("System: p + 40Ca at E_lab = $E_lab MeV (E_cm = $E_cm MeV)")
    println()

    # Create potential with Coulomb
    pot = OpticalPotential(
        V_v=53.3, r_v=1.25, a_v=0.65,
        W_v=0.0, r_wv=1.25, a_wv=0.65,
        V_s=0.0, r_s=1.25, a_s=0.65,
        W_s=13.5, r_ws=1.25, a_ws=0.47,
        V_so=7.5, r_so=1.10, a_so=0.65,
        W_so=0.0, r_wso=1.10, a_wso=0.65,
        r_c=1.25,
        Z_proj=Z_proj, A_proj=A_proj,
        Z_targ=Z_targ, A_targ=A_targ,
        A1=0.0, A2=0.0  # COLOSS convention
    )

    # Physical parameters
    HBARC = 197.327
    amu = 931.5
    m_proj = A_proj * amu
    m_targ = A_targ * amu
    μ = m_proj * m_targ / (m_proj + m_targ)
    k = sqrt(2 * μ * E_cm) / HBARC
    η = SLAM.coulomb_eta(Z_proj, Z_targ, μ, E_cm)
    e2 = 1.44  # MeV·fm

    println("Physical parameters:")
    println("  μ = $(round(μ, digits=2)) MeV/c²")
    println("  k = $(round(k, digits=6)) fm⁻¹")
    println("  η = $(round(η, digits=4)) (Sommerfeld parameter)")
    println()

    # Test short-range potential at several radii
    println("Short-range potential test:")
    println("-" ^ 70)
    println("  r (fm)   V_nuc (MeV)    V_coul_finite   V_coul_point   V_short")
    println("-" ^ 70)

    a13 = A_targ^(1/3)  # COLOSS convention
    R_c = 1.25 * a13
    println("  Coulomb radius R_c = $(round(R_c, digits=3)) fm")
    println()

    for r in [1.0, 3.0, 5.0, 8.0, 12.0]
        V_nuc = SLAM.evaluate_potential(pot, r, 0, 0.5)
        V_coul_finite = SLAM.evaluate_coulomb(pot, r)
        V_coul_point = e2 * Z_proj * Z_targ / r
        V_short = SLAM.evaluate_short_range(pot, r, 0, 0.5)

        # Expected: V_short = V_nuc + V_coul_finite - V_coul_point
        V_short_expected = V_nuc + V_coul_finite - V_coul_point

        @printf("  %4.1f  %12.4f+%.4fi  %10.4f  %10.4f  %12.4f+%.4fi\n",
                r, real(V_nuc), abs(imag(V_nuc)),
                V_coul_finite, V_coul_point,
                real(V_short), abs(imag(V_short)))

        # Verify the calculation
        if abs(V_short - V_short_expected) > 1e-10
            println("    ERROR: V_short mismatch!")
        end
    end
    println()

    # Note about the short-range potential
    println("Physical interpretation:")
    println("  - V_coul_finite = e²Z₁Z₂/r for r ≥ R_c (point Coulomb)")
    println("  - V_coul_finite = e²Z₁Z₂/(2R_c)*(3-(r/R_c)²) for r < R_c")
    println("  - V_short = V_nuc + (V_coul_finite - V_coul_point)")
    println("  - For r ≥ R_c: V_short = V_nuc (Coulomb cancels)")
    println("  - For r < R_c: V_short = V_nuc + Coulomb_correction")
    println()

    # Now test the scattering solver
    println("=" ^ 70)
    println("SCATTERING CALCULATION (l=0)")
    println("=" ^ 70)

    l = 0
    j = 0.5
    N = 60
    R = 20.0

    prob = ScatteringProblem(pot, E_cm, l; j=j, N=N, R=R)
    result = solve_scattering(prob)

    println()
    println("S-matrix: $(result.S_matrix)")
    println("|S| = $(abs(result.S_matrix))")
    println()

    # For a physical S-matrix with absorption, |S| should be < 1
    if abs(result.S_matrix) > 1.01
        println("WARNING: |S| > 1 indicates numerical issues")
        println("This may be due to incomplete Coulomb treatment at the boundary")
    else
        println("S-matrix magnitude is physical (|S| ≤ 1)")
    end

    println()
    println("Note: The current Method 5 with theta=0 (no complex scaling)")
    println("may have issues with Coulomb scattering near the boundary.")
    println("For charged particle scattering, complex scaling (Method 1/2)")
    println("or larger integration radius may be needed.")
end

main()
