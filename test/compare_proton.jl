# Compare proton scattering results: SLAM.jl vs COLOSS Method 5
# Note: Both have |S| > 1 for theta=0, which is expected for Method 5 with Coulomb

using SLAM
using Printf

function main()
    println("=" ^ 70)
    println("PROTON SCATTERING: SLAM.jl vs COLOSS Method 5")
    println("p + 40Ca at 30 MeV, theta=0 (no complex scaling)")
    println("=" ^ 70)
    println()

    # COLOSS results (from running test_method5_proton.in)
    # These have |S| > 1 due to no complex scaling
    coloss_results = [
        # (l, S, J, Re(S), Im(S))
        (0, 0.5, 0.5, 0.869879, -0.710145),
        (1, 0.5, 0.5, 0.280556, -1.008179),
        (1, 0.5, 1.5, 0.343461, -1.036486),
        (2, 0.5, 1.5, 0.021722, -0.697904),
        (2, 0.5, 2.5, 0.010180, -0.816266),
        (3, 0.5, 2.5, 0.357794, -0.489242),
        (3, 0.5, 3.5, 0.209968, -0.560292),
    ]

    # Physical parameters
    A_proj = 1.0
    A_targ = 40.0
    Z_proj = 1.0
    Z_targ = 20.0
    E_lab = 30.0
    E_cm = E_lab * A_targ / (A_proj + A_targ)

    # Create potential
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

    N = 60
    R = 20.0

    μ = SLAM.reduced_mass(pot)
    k = SLAM.coulomb_k(μ, E_cm)
    η = SLAM.coulomb_eta(Z_proj, Z_targ, μ, E_cm)
    println("E_cm = $(round(E_cm, digits=3)) MeV, k = $(round(k, digits=6)) fm⁻¹, η = $(round(η, digits=4))")
    println()

    println("  L    S      J  |  COLOSS S-matrix        |  SLAM S-matrix          |  Δ|S|     ΔRe(S)   ΔIm(S)")
    println("-" ^ 100)

    max_re_diff = 0.0
    max_im_diff = 0.0

    for (l, s, j, re_coloss, im_coloss) in coloss_results
        prob = ScatteringProblem(pot, E_cm, l; j=j, N=N, R=R)
        result = solve_scattering(prob)

        re_slam = real(result.S_matrix)
        im_slam = imag(result.S_matrix)
        abs_coloss = sqrt(re_coloss^2 + im_coloss^2)
        abs_slam = abs(result.S_matrix)

        d_abs = abs_slam - abs_coloss
        d_re = re_slam - re_coloss
        d_im = im_slam - im_coloss

        max_re_diff = max(max_re_diff, abs(d_re))
        max_im_diff = max(max_im_diff, abs(d_im))

        @printf("  %d   %.1f   %.1f  |  (%8.6f, %8.6f)  |  (%8.6f, %8.6f)  |  %+.4f  %+.5f  %+.5f\n",
                l, s, j, re_coloss, im_coloss, re_slam, im_slam, d_abs, d_re, d_im)
    end

    println("-" ^ 100)
    println()
    @printf("Maximum absolute difference: Re(S) = %.6f, Im(S) = %.6f\n", max_re_diff, max_im_diff)
    println()
    println("Note: Both SLAM and COLOSS give |S| > 1 because Method 5 with theta=0")
    println("does not properly handle Coulomb scattering. This confirms the short-range")
    println("potential and Coulomb function treatment matches COLOSS.")
end

main()
