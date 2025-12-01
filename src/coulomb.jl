"""
    Coulomb wave functions for scattering calculations.

Uses the COUL90 Fortran library (A.R. Barnett) via ccall for accurate
evaluation of regular (F) and irregular (G) Coulomb functions,
as well as the outgoing Hankel function H⁺ = G + iF.

Reference:
- A.R. Barnett, CPC 21 (1981) 297-314
- NIST DLMF: https://dlmf.nist.gov/33
"""

# Path to the compiled Fortran library
const _libcoul90_path = Ref{String}("")
const _libcoul90_loaded = Ref{Bool}(false)

function __init_coulomb__()
    # Try to find and load the library
    libname = Sys.isapple() ? "libcoul90.dylib" : Sys.iswindows() ? "libcoul90.dll" : "libcoul90.so"
    libpath = joinpath(@__DIR__, "..", "deps", libname)

    if isfile(libpath)
        _libcoul90_path[] = libpath
        _libcoul90_loaded[] = true
    else
        @warn "Coulomb library not found at $libpath. Using pure Julia fallback. Run the build script to compile."
        _libcoul90_loaded[] = false
    end
end

"""
    coulomb_eta(Z1::Float64, Z2::Float64, μ::Float64, E::Float64) -> Float64

Compute the Sommerfeld parameter η for Coulomb scattering.

η = Z₁ Z₂ e² μ / (ℏ² k)

where k = √(2μE)/ℏ is the wave number.

# Arguments
- `Z1, Z2`: Charge numbers of the two particles
- `μ`: Reduced mass in MeV/c²
- `E`: Center of mass energy in MeV

# Returns
- `Float64`: Sommerfeld parameter η
"""
function coulomb_eta(Z1::Float64, Z2::Float64, μ::Float64, E::Float64)
    α = 1.0 / 137.036  # Fine structure constant
    ℏc = 197.327       # MeV·fm
    k = sqrt(2 * μ * E) / ℏc
    return Z1 * Z2 * α * μ / (ℏc * k)
end

"""
    coulomb_k(μ::Float64, E::Float64) -> Float64

Compute the wave number k for given reduced mass and energy.

k = √(2μE) / ℏc  [in fm⁻¹]
"""
function coulomb_k(μ::Float64, E::Float64)
    ℏc = 197.327  # MeV·fm
    return sqrt(2 * μ * E) / ℏc
end

"""
    coulomb_FG(l::Int, η::Float64, ρ::Float64) -> Tuple{Float64, Float64, Float64, Float64}

Compute Coulomb functions F_l, G_l and their derivatives F'_l, G'_l.

If the Fortran COUL90 library is available, uses that for high accuracy.
Otherwise falls back to a pure Julia implementation.

# Arguments
- `l`: Angular momentum quantum number (≥ 0)
- `η`: Sommerfeld parameter
- `ρ`: Dimensionless radius ρ = kr (> 0)

# Returns
- `(F, G, Fp, Gp)`: F_l, G_l, F'_l, G'_l at ρ
"""
function coulomb_FG(l::Int, η::Float64, ρ::Float64)
    if _libcoul90_loaded[]
        return _coulomb_FG_fortran(l, η, ρ)
    else
        return _coulomb_FG_julia(l, η, ρ)
    end
end

function _coulomb_FG_fortran(l::Int, η::Float64, ρ::Float64)
    fc = Ref{Float64}(0.0)
    gc = Ref{Float64}(0.0)
    fcp = Ref{Float64}(0.0)
    gcp = Ref{Float64}(0.0)
    ifail = Ref{Cint}(0)

    ccall((:coul90_single, _libcoul90_path[]), Cvoid,
          (Float64, Float64, Cint, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Cint}),
          ρ, η, l, fc, gc, fcp, gcp, ifail)

    if ifail[] != 0
        @warn "COUL90 returned error code $(ifail[])"
    end

    return (fc[], gc[], fcp[], gcp[])
end

# Pure Julia fallback implementation using continued fractions
function _coulomb_FG_julia(l::Int, η::Float64, ρ::Float64)
    if ρ ≤ 0
        error("ρ must be positive")
    end

    lf = Float64(l)

    # Compute f = F'/F using continued fraction CF1
    f = _cf1_coulomb(l, η, ρ)

    # Compute p + iq = H'/H using continued fraction CF2
    p, q = _cf2_coulomb(l, η, ρ)

    w = f - p
    den = w^2 + q^2

    F = 1.0 / sqrt(den)
    if ρ < 0.1 * (lf + 1)
        F = abs(F)
    end

    G = w * F / q
    Fp = f * F
    Gp = p * G + q * F

    return (F, G, Fp, Gp)
end

function _cf1_coulomb(l::Int, η::Float64, ρ::Float64)
    lf = Float64(l)
    small = 1e-50
    acc = 1e-15

    f0 = (lf + 1) / ρ + η / (lf + 1)
    if abs(f0) < small
        f0 = small
    end

    f = f0
    C = f0
    D = 0.0

    for n in 1:10000
        nf = Float64(n)
        lpn = lf + nf

        an = -(η^2 + lpn^2) * ((lpn)^2 - lf^2) / ((2*lpn - 1) * (2*lpn + 1))
        bn = 2 * (lpn / ρ + η / lpn)

        D = bn + an * D
        if abs(D) < small
            D = small
        end
        D = 1.0 / D

        C = bn + an / C
        if abs(C) < small
            C = small
        end

        delta = C * D
        f = f * delta

        if abs(delta - 1.0) < acc
            return f
        end
    end

    @warn "CF1 did not converge"
    return f
end

function _cf2_coulomb(l::Int, η::Float64, ρ::Float64)
    lf = Float64(l)
    small = 1e-50
    acc = 1e-15

    p0 = 0.0
    q0 = 1.0 - η / ρ

    ar = 2 * (η - lf)
    ai = 2.0
    br = 2 * ρ
    bi = 2 * η

    denom = br^2 + bi^2
    dr = br / denom
    di = -bi / denom

    cr = br
    ci = bi

    pr = p0 + ar * dr - ai * di
    pi_val = q0 + ar * di + ai * dr

    for n in 1:10000
        nf = Float64(n)

        ar = -(lf + nf) * (lf - nf + 1) + η^2
        ai = η * (2 * nf - 1)
        br = 2 * (ρ - η + nf)
        bi = 2 * η

        denom_r = br + ar * dr - ai * di
        denom_i = bi + ar * di + ai * dr
        denom2 = denom_r^2 + denom_i^2
        if denom2 < small^2
            denom2 = small^2
        end
        dr_new = denom_r / denom2
        di_new = -denom_i / denom2

        cdenom = cr^2 + ci^2
        if cdenom < small^2
            cdenom = small^2
        end
        acr = ar * cr + ai * ci
        aci = ai * cr - ar * ci
        cr_new = br + acr / cdenom
        ci_new = bi + aci / cdenom

        delta_r = cr_new * dr_new - ci_new * di_new
        delta_i = cr_new * di_new + ci_new * dr_new

        pr_new = pr * delta_r - pi_val * delta_i
        pi_new = pr * delta_i + pi_val * delta_r

        if abs(delta_r - 1.0) + abs(delta_i) < acc
            return (pr_new, pi_new)
        end

        pr = pr_new
        pi_val = pi_new
        dr = dr_new
        di = di_new
        cr = cr_new
        ci = ci_new
    end

    @warn "CF2 did not converge"
    return (pr, pi_val)
end

"""
    coulomb_F(l::Int, η::Float64, ρ::Float64) -> Float64

Compute the regular Coulomb function F_l(η, ρ).
"""
function coulomb_F(l::Int, η::Float64, ρ::Float64)
    F, G, Fp, Gp = coulomb_FG(l, η, ρ)
    return F
end

"""
    coulomb_G(l::Int, η::Float64, ρ::Float64) -> Float64

Compute the irregular Coulomb function G_l(η, ρ).
"""
function coulomb_G(l::Int, η::Float64, ρ::Float64)
    F, G, Fp, Gp = coulomb_FG(l, η, ρ)
    return G
end

"""
    coulomb_H_plus(l::Int, η::Float64, ρ::Float64) -> ComplexF64

Compute the outgoing Coulomb-Hankel function H⁺_l = G_l + i*F_l.
"""
function coulomb_H_plus(l::Int, η::Float64, ρ::Float64)
    F, G, Fp, Gp = coulomb_FG(l, η, ρ)
    return complex(G, F)
end

"""
    coulomb_H_plus_derivative(l::Int, η::Float64, ρ::Float64) -> ComplexF64

Compute the derivative of the outgoing Coulomb-Hankel function.
"""
function coulomb_H_plus_derivative(l::Int, η::Float64, ρ::Float64)
    F, G, Fp, Gp = coulomb_FG(l, η, ρ)
    return complex(Gp, Fp)
end

"""
    coulomb_gamma_s(l::Int, η::Float64, ρ::Float64) -> ComplexF64

Compute the logarithmic derivative of the outgoing wave at r = R:
γ_s = d/dr[H⁺_l(kr)] / H⁺_l(kr) = k * H'⁺_l(ρ) / H⁺_l(ρ)

where H'⁺_l is the derivative with respect to ρ = kr.

This is used for the outgoing wave boundary condition at r = R.

Note: The k factor is included to convert from dρ to dr derivative.
"""
function coulomb_gamma_s(l::Int, η::Float64, ρ::Float64, k::Float64)
    F, G, Fp, Gp = coulomb_FG(l, η, ρ)
    H_plus = complex(G, F)
    H_plus_prime_ρ = complex(Gp, Fp)  # derivative w.r.t. ρ
    # Convert to derivative w.r.t. r: d/dr = k * d/dρ
    return k * H_plus_prime_ρ / H_plus
end

# Keep old signature for compatibility (without k factor, for ρ-derivative)
function coulomb_gamma_s(l::Int, η::Float64, ρ::Float64)
    F, G, Fp, Gp = coulomb_FG(l, η, ρ)
    H_plus = complex(G, F)
    H_plus_prime = complex(Gp, Fp)
    return H_plus_prime / H_plus
end

"""
    coulomb_sigma(l::Int, η::Float64) -> Float64

Compute the Coulomb phase shift σ_l = arg(Γ(l+1+iη)).
"""
function coulomb_sigma(l::Int, η::Float64)
    z = complex(Float64(l + 1), η)
    return imag(loggamma(z))
end
