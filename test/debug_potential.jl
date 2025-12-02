# Debug: Detailed potential comparison
using SLAM

function main()
    println("=" ^ 70)
    println("Potential value debug")
    println("=" ^ 70)

    # Test parameters
    pot = OpticalPotential(
        V_v=46.553, r_v=1.185, a_v=0.672,
        W_v=1.777, r_wv=1.185, a_wv=0.672,
        V_s=0.0, r_s=1.288, a_s=0.538,
        W_s=7.182, r_ws=1.288, a_ws=0.538,
        V_so=5.343, r_so=0.996, a_so=0.590,
        W_so=-0.110, r_wso=0.996, a_wso=0.590,
        r_c=0.0,
        Z_proj=0.0, A_proj=1.0,
        Z_targ=0.0, A_targ=40.0
    )

    r = 5.0
    l = 0
    j = 0.5

    a13 = pot.A_proj^(1/3) + pot.A_targ^(1/3)
    println("a13 = $(a13)")

    R_v = pot.r_v * a13
    R_wv = pot.r_wv * a13
    R_ws = pot.r_ws * a13
    R_so = pot.r_so * a13
    R_wso = pot.r_wso * a13

    println("R_v = $(R_v)")
    println("R_wv = $(R_wv)")
    println("R_ws = $(R_ws)")
    println("R_so = $(R_so)")
    println("R_wso = $(R_wso)")
    println()

    # Woods-Saxon form factors
    f_v = woods_saxon(r, R_v, pot.a_v)
    f_wv = woods_saxon(r, R_wv, pot.a_wv)
    println("f_v(r=$r) = $(f_v)")
    println("f_wv(r=$r) = $(f_wv)")
    println()

    # Volume term
    V_volume = -pot.V_v * f_v - im * pot.W_v * f_wv
    println("V_volume = -$(pot.V_v) * $(f_v) - i*$(pot.W_v) * $(f_wv)")
    println("         = $(V_volume)")
    println()

    # Surface term
    df_ws = woods_saxon_derivative(r, R_ws, pot.a_ws)
    println("df_ws(r=$r) = $(df_ws)")
    V_surface = im * 4.0 * pot.W_s * pot.a_ws * df_ws
    println("V_surface = 4i * $(pot.W_s) * $(pot.a_ws) * $(df_ws)")
    println("          = $(V_surface)")
    println()

    # Spin-orbit
    s = 0.5
    ls = 0.5 * (j*(j+1) - l*(l+1) - s*(s+1))
    println("ls = 0.5 * ($(j)*($(j)+1) - $(l)*($(l)+1) - 0.5*1.5) = $(ls)")

    df_so = woods_saxon_derivative(r, R_so, pot.a_so) / r
    df_wso = woods_saxon_derivative(r, R_wso, pot.a_wso) / r
    println("df_so/r = $(df_so)")
    println("df_wso/r = $(df_wso)")

    V_spinorbit = 2.0 * pot.V_so * df_so * ls + im * 2.0 * pot.W_so * df_wso * ls
    println("V_spinorbit = 2 * $(pot.V_so) * $(df_so) * $(ls) + 2i * $(pot.W_so) * $(df_wso) * $(ls)")
    println("            = $(V_spinorbit)")
    println()

    # Total
    V_total = V_volume + V_surface + V_spinorbit
    println("V_total = V_volume + V_surface + V_spinorbit")
    println("        = $(V_total)")
    println()

    # Now using evaluate_potential
    V_from_func = evaluate_potential(pot, r, l, j)
    println("evaluate_potential(pot, r=$r, l=$l, j=$j) = $(V_from_func)")
    println()

    # Short-range
    V_short = evaluate_short_range(pot, r, l, j)
    println("evaluate_short_range(pot, r=$r, l=$l, j=$j) = $(V_short)")
    println()

    # Convert to U
    HBARC = 197.327
    A_proj = 1.0
    A_targ = 40.0
    amu = 931.5
    m_proj = A_proj * amu
    m_targ = A_targ * amu
    μ = m_proj * m_targ / (m_proj + m_targ)
    coeff_kin = HBARC^2 / (2.0 * μ)

    println("μ = $(μ) MeV/c²")
    println("coeff_kin = ℏ²/(2μ) = $(coeff_kin) MeV·fm²")

    U = V_short / coeff_kin
    println()
    println("U = V_short / coeff_kin = $(U) fm⁻²")

    println()
    println("=" ^ 70)
    println("Potential calculation details:")
    println("=" ^ 70)
    println()

    # Verify Woods-Saxon derivative formula
    V_central = -pot.V_v * f_v - im * pot.W_v * f_wv
    println("V_central = $V_central")

    # Woods-Saxon derivative: df/dr = -exp((r-R)/a)/a/(1+exp((r-R)/a))^2
    # This is equivalent to: df/dr = -f*(1-f)/a where f = 1/(1+exp(x))
    x_ws = (r - R_ws) / pot.a_ws
    exp_x = exp(x_ws)
    df_ws_formula = -exp_x / pot.a_ws / (1 + exp_x)^2

    println("df_ws (woods_saxon_derivative) = $(df_ws)")
    println("df_ws (formula) = $(df_ws_formula)")

    V_surface_calc = 4.0 * im * pot.W_s * pot.a_ws * df_ws_formula
    println("V_surface = $V_surface_calc")

    # Spin-orbit
    x_so = (r - R_so) / pot.a_so
    exp_so = exp(x_so)
    df_so_formula = -exp_so / pot.a_so / (1 + exp_so)^2 / r

    x_wso = (r - R_wso) / pot.a_wso
    exp_wso = exp(x_wso)
    df_wso_formula = -exp_wso / pot.a_wso / (1 + exp_wso)^2 / r

    println("df_so (SLAM) = $(df_so)")
    println("df_so (formula) = $(df_so_formula)")

    V_so_calc = 2.0 * pot.V_so * df_so_formula * ls + 2.0 * im * pot.W_so * df_wso_formula * ls
    println("V_so = $V_so_calc")

    V_nuc_calc = V_central + V_surface_calc + V_so_calc
    println()
    println("V_nuc (calculated) = $V_nuc_calc")
    println("V_nuc (SLAM)       = $V_from_func")
    println("Difference = $(V_nuc_calc - V_from_func)")
end

main()
