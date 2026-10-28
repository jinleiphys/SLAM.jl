# Neutron scattering validation test for SLAM.jl
# Test case: n + 40Target at Elab = 20 MeV (no Coulomb)

using SLAM
using Printf

function main()
    println("=" ^ 70)
    println("SLAM.jl Neutron Scattering Validation")
    println("Test: n + 40Target at Elab = 20 MeV (no Coulomb)")
    println("=" ^ 70)

    # Parameters
    # zp=0, zt=0 (no Coulomb)
    # massp=1, masst=40
    # elab=20, sp=0.5
    # nr=60, Rmax=20

    # Convert Elab to Ecm
    E_lab = 20.0
    A_proj = 1.0
    A_targ = 40.0
    E_cm = E_lab * A_targ / (A_proj + A_targ)

    println("Lab Energy: $E_lab MeV")
    println("CM Energy: $E_cm MeV")
    println()

    # Create optical potential (no Coulomb: zp=0, zt=0)
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

    # Reference results: COLOSS (src/rot_potential.f, l·σ spin-orbit), run 2026-10-28 on
    # COLOSS/test/conv_m1_noCoul_ref.in (nr=60, Rmax=50, ctheta=6). Expect agreement to ~3e-3.
    # The previous table was not reproducible by COLOSS and matched the old half-strength spin-orbit.
    reference_results = [
        # (l, s, j, Re(S), Im(S), σ_partial)
        (0, 0.5, 0.5, 0.349785, 0.223248, 28.5536),
        (1, 0.5, 0.5, 0.508540, 0.041798, 25.5123),
        (1, 0.5, 1.5, 0.479625, 0.195694, 50.4745),
        (2, 0.5, 1.5, 0.321685, -0.152676, 60.2391),
        (2, 0.5, 2.5, 0.392000, 0.086104, 86.8107),
        (3, 0.5, 2.5, -0.032093, -0.496280, 77.8861),
        (3, 0.5, 3.5, 0.317072, -0.310444, 110.8037),
        (4, 0.5, 3.5, -0.134345, -0.189786, 130.5120),
        (4, 0.5, 4.5, 0.172014, -0.508817, 122.7115),
        (5, 0.5, 4.5, 0.448928, 0.217776, 129.5274),
        (5, 0.5, 5.5, 0.219609, 0.066820, 196.0524),
        (6, 0.5, 5.5, 0.878904, 0.126747, 43.7637),
        (6, 0.5, 6.5, 0.839524, 0.167095, 64.5348),
        (7, 0.5, 6.5, 0.977510, 0.037757, 10.3943),
        (7, 0.5, 7.5, 0.975108, 0.047042, 12.9562),
        (8, 0.5, 7.5, 0.995748, 0.009952, 2.3144),
        (8, 0.5, 8.5, 0.995614, 0.011765, 2.6742),
        (9, 0.5, 8.5, 0.999183, 0.002546, 0.5053),
        (9, 0.5, 9.5, 0.999178, 0.002913, 0.5636),
        (10, 0.5, 9.5, 0.999841, 0.000644, 0.1093),
    ]

    println("  L    S      J  |  Reference S-matrix     |  SLAM S-matrix          |  Δ|S|     ΔRe(S)   ΔIm(S)")
    println("-" ^ 100)

    total_diff_re = 0.0
    total_diff_im = 0.0
    n_channels = 0

    for (l, s, j, re_ref, im_ref, σ_ref) in reference_results
        # Solve with SLAM
        prob = ScatteringProblem(pot, E_cm, l; j=j, N=60, R=20.0)
        result = solve_scattering(prob)

        S_slam = result.S_matrix
        re_slam = real(S_slam)
        im_slam = imag(S_slam)

        S_ref = complex(re_ref, im_ref)

        # Compute differences
        diff_re = re_slam - re_ref
        diff_im = im_slam - im_ref
        diff_abs = abs(S_slam) - abs(S_ref)

        total_diff_re += abs(diff_re)
        total_diff_im += abs(diff_im)
        n_channels += 1

        @printf("  %d   %.1f   %.1f  |  (%8.5f, %8.5f)  |  (%8.5f, %8.5f)  |  %+.4f  %+.5f  %+.5f\n",
                l, s, j, re_ref, im_ref, re_slam, im_slam, diff_abs, diff_re, diff_im)
    end

    println("-" ^ 100)
    println()
    @printf("Mean absolute difference: Re(S) = %.6f, Im(S) = %.6f\n",
            total_diff_re / n_channels, total_diff_im / n_channels)
    println()

    # Print a few detailed comparisons
    println("\nDetailed comparison for l=0:")
    prob = ScatteringProblem(pot, E_cm, 0; j=0.5, N=60, R=20.0)
    result = solve_scattering(prob)
    println("  SLAM:      S = $(result.S_matrix)")
    println("  Reference: S = $(complex(0.349464, 0.223937))")
    @printf("  |S|_SLAM = %.6f, |S|_ref = %.6f\n",
            abs(result.S_matrix), abs(complex(0.349464, 0.223937)))
end

main()
