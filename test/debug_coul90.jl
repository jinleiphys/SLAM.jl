# Debug: Compare COUL90 Fortran output with Julia fallback
using SLAM

function main()
    println("=" ^ 70)
    println("COUL90 comparison: Fortran vs Julia fallback")
    println("=" ^ 70)

    # Check if Fortran library is loaded
    println("Library loaded: $(SLAM._libcoul90_loaded[])")
    println("Library path: $(SLAM._libcoul90_path[])")
    println()

    # Test parameters
    l = 0
    η = 0.0
    ρ = 19.087148  # k*R

    # Force Fortran
    F_fort, G_fort, Fp_fort, Gp_fort = SLAM._coulomb_FG_fortran(l, η, ρ)

    # Force Julia fallback
    F_julia, G_julia, Fp_julia, Gp_julia = SLAM._coulomb_FG_julia(l, η, ρ)

    println("l=$l, η=$η, ρ=$ρ:")
    println()
    println("  F (Fortran): $F_fort")
    println("  F (Julia):   $F_julia")
    println("  F (exact):   $(sin(ρ))")
    println()
    println("  G (Fortran): $G_fort")
    println("  G (Julia):   $G_julia")
    println("  G (exact):   $(cos(ρ))")
    println()
    println("  F' (Fortran): $Fp_fort")
    println("  F' (Julia):   $Fp_julia")
    println("  F' (exact):   $(cos(ρ))")
    println()
    println("  G' (Fortran): $Gp_fort")
    println("  G' (Julia):   $Gp_julia")
    println("  G' (exact):   $(-sin(ρ))")
    println()

    # Check H+ and gamma_s
    H_plus_fort = complex(G_fort, F_fort)
    H_plus_julia = complex(G_julia, F_julia)
    H_plus_exact = exp(im * ρ)

    println("H⁺ = G + iF:")
    println("  Fortran: $H_plus_fort")
    println("  Julia:   $H_plus_julia")
    println("  Exact:   $H_plus_exact")
    println()

    # Now test what coulomb_FG returns (uses either Fortran or Julia depending on availability)
    F, G, Fp, Gp = coulomb_FG(l, η, ρ)
    println("coulomb_FG returns:")
    println("  F=$F, G=$G, Fp=$Fp, Gp=$Gp")
end

main()
