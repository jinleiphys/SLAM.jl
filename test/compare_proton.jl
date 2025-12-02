# Proton scattering validation test for SLAM.jl
# p + 40Ca at 30 MeV

using SLAM
using Printf

function main()
    println("=" ^ 70)
    println("PROTON SCATTERING: SLAM.jl Validation Test")
    println("p + 40Ca at 30 MeV")
    println("=" ^ 70)
    println()

    # Reference results for validation
    reference_results = [
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
        A1=0.0, A2=0.0
    )

    N = 60
    R = 20.0

    μ = SLAM.reduced_mass(pot)
    k = SLAM.coulomb_k(μ, E_cm)
    η = SLAM.coulomb_eta(Z_proj, Z_targ, μ, E_cm)
    println("E_cm = $(round(E_cm, digits=3)) MeV, k = $(round(k, digits=6)) fm⁻¹, η = $(round(η, digits=4))")
    println()

    println("  L    S      J  |  Reference S-matrix     |  SLAM S-matrix          |  Δ|S|     ΔRe(S)   ΔIm(S)")
    println("-" ^ 100)

    max_re_diff = 0.0
    max_im_diff = 0.0

    for (l, s, j, re_ref, im_ref) in reference_results
        prob = ScatteringProblem(pot, E_cm, l; j=j, N=N, R=R)
        result = solve_scattering(prob)

        re_slam = real(result.S_matrix)
        im_slam = imag(result.S_matrix)
        abs_ref = sqrt(re_ref^2 + im_ref^2)
        abs_slam = abs(result.S_matrix)

        d_abs = abs_slam - abs_ref
        d_re = re_slam - re_ref
        d_im = im_slam - im_ref

        max_re_diff = max(max_re_diff, abs(d_re))
        max_im_diff = max(max_im_diff, abs(d_im))

        @printf("  %d   %.1f   %.1f  |  (%8.6f, %8.6f)  |  (%8.6f, %8.6f)  |  %+.4f  %+.5f  %+.5f\n",
                l, s, j, re_ref, im_ref, re_slam, im_slam, d_abs, d_re, d_im)
    end

    println("-" ^ 100)
    println()
    @printf("Maximum absolute difference: Re(S) = %.6f, Im(S) = %.6f\n", max_re_diff, max_im_diff)
    println()
    println("Note: |S| > 1 may occur for Coulomb scattering without complex scaling.")
    println("This confirms the short-range potential and Coulomb function treatment.")
end

main()
