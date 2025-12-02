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

    # Reference results for validation
    reference_results = [
        # (l, s, j, Re(S), Im(S), σ_partial)
        (0, 0.5, 0.5, 0.349464, 0.223937, 28.5507),
        (1, 0.5, 0.5, 0.507241, 0.094951, 25.3072),
        (1, 0.5, 1.5, 0.490236, 0.171754, 50.3714),
        (2, 0.5, 1.5, 0.356919, -0.085196, 59.6970),
        (2, 0.5, 2.5, 0.389119, 0.036551, 87.6725),
        (3, 0.5, 2.5, 0.090488, -0.482903, 78.5007),
        (3, 0.5, 3.5, 0.264857, -0.377562, 108.6249),
        (4, 0.5, 3.5, -0.105955, -0.313711, 122.8444),
        (4, 0.5, 4.5, 0.065042, -0.482411, 131.5992),
        (5, 0.5, 4.5, 0.387865, 0.206886, 139.1375),
        (5, 0.5, 5.5, 0.266528, 0.133520, 188.5664),
        (6, 0.5, 5.5, 0.870474, 0.137035, 46.2543),
        (6, 0.5, 6.5, 0.850974, 0.157380, 60.6222),
        (7, 0.5, 6.5, 0.976940, 0.040114, 10.6187),
        (7, 0.5, 7.5, 0.975742, 0.044760, 12.6725),
        (8, 0.5, 7.5, 0.995739, 0.010433, 2.3167),
        (8, 0.5, 8.5, 0.995671, 0.011338, 2.6419),
        (9, 0.5, 8.5, 0.999187, 0.002644, 0.5027),
        (9, 0.5, 9.5, 0.999185, 0.002828, 0.5596),
        (10, 0.5, 9.5, 0.999842, 0.000664, 0.1092),
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
