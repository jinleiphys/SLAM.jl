# Debug: Verify Coulomb functions for eta=0 case
using SLAM

function main()
    println("=" ^ 60)
    println("Coulomb function verification (η=0 case)")
    println("=" ^ 60)
    println()

    η = 0.0

    # For η=0, l=0:
    # F_0(ρ) = sin(ρ)
    # G_0(ρ) = cos(ρ)
    # F'_0(ρ) = cos(ρ)
    # G'_0(ρ) = -sin(ρ)

    println("l=0, η=0:")
    println("ρ        F_0 (SLAM)    sin(ρ)        G_0 (SLAM)    cos(ρ)      ")
    println("-" ^ 70)

    for ρ in [0.5, 1.0, 5.0, 10.0, 19.0869]
        F, G, Fp, Gp = coulomb_FG(0, η, ρ)

        exact_F = sin(ρ)
        exact_G = cos(ρ)
        exact_Fp = cos(ρ)
        exact_Gp = -sin(ρ)

        println("ρ=$ρ:")
        println("  F = $F, exact = $exact_F, diff = $(F - exact_F)")
        println("  G = $G, exact = $exact_G, diff = $(G - exact_G)")
        println("  F' = $Fp, exact = $exact_Fp, diff = $(Fp - exact_Fp)")
        println("  G' = $Gp, exact = $exact_Gp, diff = $(Gp - exact_Gp)")
        println()
    end

    # Test H+ and gamma_s
    println("=" ^ 60)
    println("H+ and γ_s verification")
    println("=" ^ 60)

    k = 0.954357  # from debug output
    R = 20.0
    ρ_R = k * R  # ≈ 19.0869

    F, G, Fp, Gp = coulomb_FG(0, η, ρ_R)

    H_plus = complex(G, F)
    H_plus_prime_ρ = complex(Gp, Fp)  # w.r.t. ρ
    H_plus_prime_r = k * H_plus_prime_ρ  # w.r.t. r

    γ_s_ρ = H_plus_prime_ρ / H_plus  # derivative w.r.t. ρ
    γ_s_r = H_plus_prime_r / H_plus  # derivative w.r.t. r

    println("At ρ = kR = $ρ_R:")
    println("  H⁺ = G + iF = $H_plus")
    println("  H⁺' (w.r.t. ρ) = G' + iF' = $H_plus_prime_ρ")
    println("  H⁺' (w.r.t. r) = k*(G' + iF') = $H_plus_prime_r")
    println()
    println("  γ_s (w.r.t. ρ) = H⁺'/H⁺ = $γ_s_ρ")
    println("  γ_s (w.r.t. r) = k*H⁺'/H⁺ = $γ_s_r")
    println()

    # Expected for η=0, l=0:
    # H⁺ = cos(ρ) + i*sin(ρ) = exp(iρ)
    # H⁺' (w.r.t. ρ) = -sin(ρ) + i*cos(ρ) = i*exp(iρ)
    # γ_s (w.r.t. ρ) = i*exp(iρ)/exp(iρ) = i
    # γ_s (w.r.t. r) = k*i

    println("Expected:")
    println("  H⁺ = exp(iρ) = $(exp(im*ρ_R))")
    println("  γ_s (w.r.t. ρ) = i = $(im)")
    println("  γ_s (w.r.t. r) = ki = $(k*im)")

    println()
    println("Using SLAM's coulomb_gamma_s functions:")
    γ_s_slam_ρ = coulomb_gamma_s(0, η, ρ_R)
    γ_s_slam_r = coulomb_gamma_s(0, η, ρ_R, k)
    println("  γ_s (w.r.t. ρ, 3-arg) = $γ_s_slam_ρ")
    println("  γ_s (w.r.t. r, 4-arg) = $γ_s_slam_r")
end

main()
