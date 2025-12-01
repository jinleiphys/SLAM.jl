"""
    Optical model potentials for nuclear scattering calculations.

Implements Woods-Saxon form factors for volume, surface, and spin-orbit
potentials, as well as Coulomb potentials with finite-size corrections.

Reference: Nuclear Reaction Data and Theory (NRV)
"""

"""
    OpticalPotential

Structure holding optical model potential parameters.

The optical potential has the form:
U(r) = V_central(r) + V_surface(r) + V_spin_orbit(r) + V_coulomb(r)

# Fields
## Volume (Woods-Saxon)
- `V_v::Float64`: Real volume depth (MeV)
- `r_v::Float64`: Real volume radius parameter (fm)
- `a_v::Float64`: Real volume diffuseness (fm)
- `W_v::Float64`: Imaginary volume depth (MeV)
- `r_wv::Float64`: Imaginary volume radius parameter (fm)
- `a_wv::Float64`: Imaginary volume diffuseness (fm)

## Surface (derivative Woods-Saxon)
- `V_s::Float64`: Real surface depth (MeV)
- `r_s::Float64`: Real surface radius parameter (fm)
- `a_s::Float64`: Real surface diffuseness (fm)
- `W_s::Float64`: Imaginary surface depth (MeV)
- `r_ws::Float64`: Imaginary surface radius parameter (fm)
- `a_ws::Float64`: Imaginary surface diffuseness (fm)

## Spin-orbit (Thomas form)
- `V_so::Float64`: Real spin-orbit depth (MeV)
- `r_so::Float64`: Real spin-orbit radius parameter (fm)
- `a_so::Float64`: Real spin-orbit diffuseness (fm)
- `W_so::Float64`: Imaginary spin-orbit depth (MeV)
- `r_wso::Float64`: Imaginary spin-orbit radius parameter (fm)
- `a_wso::Float64`: Imaginary spin-orbit diffuseness (fm)

## Coulomb
- `r_c::Float64`: Coulomb radius parameter (fm)

## Target/Projectile info
- `Z_proj::Float64`: Projectile charge number
- `A_proj::Float64`: Projectile mass number
- `Z_targ::Float64`: Target charge number
- `A_targ::Float64`: Target mass number

## Radius calculation masses
- `A1::Float64`: Mass number for radius factor (projectile side), 0 by default
- `A2::Float64`: Mass number for radius factor (target side), defaults to A_targ

Note: Radii are calculated as R = r * (A1^(1/3) + A2^(1/3))
Following COLOSS convention: if A1=A2=0, uses A1=0 and A2=A_targ.
"""
struct OpticalPotential
    # Volume terms
    V_v::Float64
    r_v::Float64
    a_v::Float64
    W_v::Float64
    r_wv::Float64
    a_wv::Float64

    # Surface terms
    V_s::Float64
    r_s::Float64
    a_s::Float64
    W_s::Float64
    r_ws::Float64
    a_ws::Float64

    # Spin-orbit terms
    V_so::Float64
    r_so::Float64
    a_so::Float64
    W_so::Float64
    r_wso::Float64
    a_wso::Float64

    # Coulomb
    r_c::Float64

    # Particle info
    Z_proj::Float64
    A_proj::Float64
    Z_targ::Float64
    A_targ::Float64

    # Radius calculation masses (COLOSS convention)
    A1::Float64
    A2::Float64
end

"""
    OpticalPotential(; kwargs...)

Construct an optical potential with keyword arguments.
Default values are zero for all potential depths and 1.25 fm for radius parameters.

## COLOSS Convention for Radius Calculation
- If A1 and A2 are not specified (or both zero), uses A1=0 and A2=A_targ
- Radii are calculated as R = r_param * (A1^(1/3) + A2^(1/3))
"""
function OpticalPotential(;
    V_v=0.0, r_v=1.25, a_v=0.65,
    W_v=0.0, r_wv=1.25, a_wv=0.65,
    V_s=0.0, r_s=1.25, a_s=0.65,
    W_s=0.0, r_ws=1.25, a_ws=0.65,
    V_so=0.0, r_so=1.10, a_so=0.65,
    W_so=0.0, r_wso=1.10, a_wso=0.65,
    r_c=1.25,
    Z_proj=1.0, A_proj=1.0,
    Z_targ=20.0, A_targ=40.0,
    A1=nothing, A2=nothing  # Radius calculation masses
)
    # Apply COLOSS convention: if A1=A2=0 or not specified, use A1=0, A2=A_targ
    if A1 === nothing && A2 === nothing
        A1_val = 0.0
        A2_val = A_targ
    elseif A1 === nothing
        A1_val = 0.0
        A2_val = Float64(A2)
    elseif A2 === nothing
        A1_val = Float64(A1)
        A2_val = A_targ
    else
        # Both specified - apply COLOSS convention if both zero
        if A1 == 0.0 && A2 == 0.0
            A1_val = 0.0
            A2_val = A_targ
        else
            A1_val = Float64(A1)
            A2_val = Float64(A2)
        end
    end

    return OpticalPotential(
        V_v, r_v, a_v, W_v, r_wv, a_wv,
        V_s, r_s, a_s, W_s, r_ws, a_ws,
        V_so, r_so, a_so, W_so, r_wso, a_wso,
        r_c, Z_proj, A_proj, Z_targ, A_targ,
        A1_val, A2_val
    )
end

"""
    reduced_mass(pot::OpticalPotential) -> Float64

Compute the reduced mass in MeV/c² from the optical potential parameters.
"""
function reduced_mass(pot::OpticalPotential)
    # AMU to MeV/c²
    amu = 931.5  # MeV/c²
    m_proj = pot.A_proj * amu
    m_targ = pot.A_targ * amu
    return m_proj * m_targ / (m_proj + m_targ)
end

"""
    woods_saxon(r::Float64, R::Float64, a::Float64) -> Float64

Woods-Saxon form factor: f(r) = 1 / (1 + exp((r-R)/a))
"""
function woods_saxon(r::Float64, R::Float64, a::Float64)
    x = (r - R) / a
    if x > 700  # Prevent overflow
        return 0.0
    elseif x < -700
        return 1.0
    else
        return 1.0 / (1.0 + exp(x))
    end
end

"""
    woods_saxon_derivative(r::Float64, R::Float64, a::Float64) -> Float64

Derivative of Woods-Saxon form factor: df/dr = -f(1-f)/a
"""
function woods_saxon_derivative(r::Float64, R::Float64, a::Float64)
    f = woods_saxon(r, R, a)
    return -f * (1.0 - f) / a
end

"""
    evaluate_potential(pot::OpticalPotential, r::Float64, l::Int, j::Float64) -> ComplexF64

Evaluate the total optical potential at radius r for given l and j.

Returns the complex potential U(r) = V(r) - i*W(r) where V is real and W is absorptive.

# Arguments
- `pot::OpticalPotential`: Potential parameters
- `r::Float64`: Radial distance (fm)
- `l::Int`: Orbital angular momentum
- `j::Float64`: Total angular momentum

# Returns
- `ComplexF64`: Complex potential value (MeV)
"""
function evaluate_potential(pot::OpticalPotential, r::Float64, l::Int, j::Float64)
    # Compute radius factor: A1^{1/3} + A2^{1/3} (following COLOSS convention)
    a13 = pot.A1^(1/3) + pot.A2^(1/3)

    # Compute radii
    R_v = pot.r_v * a13
    R_wv = pot.r_wv * a13
    R_s = pot.r_s * a13
    R_ws = pot.r_ws * a13
    R_so = pot.r_so * a13
    R_wso = pot.r_wso * a13

    # Volume term: -V_v * f_v - i*W_v * f_wv
    f_v = woods_saxon(r, R_v, pot.a_v)
    f_wv = woods_saxon(r, R_wv, pot.a_wv)
    V_volume = -pot.V_v * f_v - im * pot.W_v * f_wv

    # Surface term: 4*a * df/dr
    # COLOSS convention: V_surface = 4*V_s*a*df + 4i*W_s*a*df
    # where df = -exp(x)/a/(1+exp(x))^2 < 0
    df_s = woods_saxon_derivative(r, R_s, pot.a_s)
    df_ws = woods_saxon_derivative(r, R_ws, pot.a_ws)
    V_surface = 4.0 * pot.V_s * pot.a_s * df_s + im * 4.0 * pot.W_s * pot.a_ws * df_ws

    # Spin-orbit term (COLOSS convention):
    # V_so = 2 * V_so * df_so/r * ls + 2i * W_so * df_wso/r * ls
    # where ls = [J(J+1) - l(l+1) - S(S+1)] / 2
    s = 0.5
    ls = 0.5 * (j*(j+1) - l*(l+1) - s*(s+1))

    V_spinorbit = 0.0 + 0.0im
    if r > 1e-6
        df_so = woods_saxon_derivative(r, R_so, pot.a_so) / r
        df_wso = woods_saxon_derivative(r, R_wso, pot.a_wso) / r
        V_spinorbit = 2.0 * pot.V_so * df_so * ls + im * 2.0 * pot.W_so * df_wso * ls
    end

    return V_volume + V_surface + V_spinorbit
end

"""
    evaluate_coulomb(pot::OpticalPotential, r::Float64) -> Float64

Evaluate the Coulomb potential at radius r.

Uses a uniformly charged sphere model for finite-size correction:
- For r ≥ r_c: V_C(r) = Z₁Z₂e²/r
- For r < r_c: V_C(r) = Z₁Z₂e²/(2r_c) * (3 - (r/r_c)²)

# Arguments
- `pot::OpticalPotential`: Potential parameters
- `r::Float64`: Radial distance (fm)

# Returns
- `Float64`: Coulomb potential (MeV)
"""
function evaluate_coulomb(pot::OpticalPotential, r::Float64)
    a13 = pot.A1^(1/3) + pot.A2^(1/3)
    R_c = pot.r_c * a13

    # e²/ℏc ≈ 1/137 and ℏc ≈ 197.3 MeV·fm
    # So e² ≈ 1.44 MeV·fm
    e2 = 1.44  # MeV·fm

    Z12 = pot.Z_proj * pot.Z_targ

    if r >= R_c
        return Z12 * e2 / r
    else
        return Z12 * e2 / (2*R_c) * (3.0 - (r/R_c)^2)
    end
end

"""
    evaluate_short_range(pot::OpticalPotential, r::Float64, l::Int, j::Float64) -> ComplexF64

Evaluate the short-range potential (nuclear + Coulomb modification) at radius r.

This is the potential that appears in the inhomogeneous term of the
scattering equation after separating out the Coulomb asymptotic behavior.

# Arguments
- `pot::OpticalPotential`: Potential parameters
- `r::Float64`: Radial distance (fm)
- `l::Int`: Orbital angular momentum
- `j::Float64`: Total angular momentum

# Returns
- `ComplexF64`: Short-range potential (MeV)
"""
function evaluate_short_range(pot::OpticalPotential, r::Float64, l::Int, j::Float64)
    # Nuclear potential
    V_nuc = evaluate_potential(pot, r, l, j)

    # For the inhomogeneous equation, we need the difference from pure Coulomb
    # In the classically forbidden region, we use the nuclear potential directly
    return V_nuc
end

"""
    evaluate_total_potential(pot::OpticalPotential, r::Float64, l::Int, j::Float64) -> ComplexF64

Evaluate the total potential (nuclear + Coulomb) at radius r.
"""
function evaluate_total_potential(pot::OpticalPotential, r::Float64, l::Int, j::Float64)
    return evaluate_potential(pot, r, l, j) + evaluate_coulomb(pot, r)
end
