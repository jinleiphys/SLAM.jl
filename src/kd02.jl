"""
    KD02 Module - Koning-Delaroche Optical Model Potential

Koning-Delaroche local and global optical model potentials for neutrons and protons.

Reference: A.J. Koning and J.P. Delaroche, Nucl. Phys. A713 (2003) 231-310
"""
module KD02

using Printf

export kd02_potential, KD02Params

"""
    KD02Params

Structure containing all optical model parameters from KD02.

# Fields
- `v`, `rv`, `av`: Real volume potential depth (MeV), radius (fm), diffuseness (fm)
- `w`, `rw`, `aw`: Imaginary volume potential depth (MeV), radius (fm), diffuseness (fm)
- `vd`, `rvd`, `avd`: Real surface potential depth (MeV), radius (fm), diffuseness (fm)
- `wd`, `rwd`, `awd`: Imaginary surface potential depth (MeV), radius (fm), diffuseness (fm)
- `vso`, `rvso`, `avso`: Real spin-orbit potential depth (MeV), radius (fm), diffuseness (fm)
- `wso`, `rwso`, `awso`: Imaginary spin-orbit potential depth (MeV), radius (fm), diffuseness (fm)
- `rc`: Coulomb radius parameter (fm)
"""
struct KD02Params
    # Real volume
    v::Float64
    rv::Float64
    av::Float64
    # Imaginary volume
    w::Float64
    rw::Float64
    aw::Float64
    # Real surface
    vd::Float64
    rvd::Float64
    avd::Float64
    # Imaginary surface
    wd::Float64
    rwd::Float64
    awd::Float64
    # Real spin-orbit
    vso::Float64
    rvso::Float64
    avso::Float64
    # Imaginary spin-orbit
    wso::Float64
    rwso::Float64
    awso::Float64
    # Coulomb
    rc::Float64
end

"""
    globalomp(k0, Z, A)

Calculate global optical model parameters.

# Arguments
- `k0::Int`: Particle type (1 = neutron, 2 = proton)
- `Z::Float64`: Target charge number
- `A::Float64`: Target mass number

# Returns
Named tuple with all geometry and energy-dependent coefficients.
"""
function globalomp(k0::Int, Z::Float64, A::Float64)
    N = A - Z

    # Common parameters for neutrons and protons
    rv = 1.3039 - 0.4054 * A^(-1/3)
    av = 0.6778 - 1.487e-4 * A
    rw = rv
    aw = av
    v4 = 7.0e-9
    w2 = 73.55 + 0.0795 * A
    rvd = 1.3424 - 0.01585 * A^(1/3)
    rwd = rvd
    vd = 0.0
    d2 = 0.0180 + 3.802e-3 / (1.0 + exp((A - 156.0) / 8.0))
    d3 = 11.5
    vso1 = 5.922 + 0.0030 * A
    vso2 = 0.0040
    rvso = 1.1854 - 0.647 * A^(-1/3)
    rwso = rvso
    avso = 0.59
    awso = avso
    wso1 = -3.1
    wso2 = 160.0

    # Particle-specific parameters
    if k0 == 1  # Neutrons
        ef = -11.2814 + 0.02646 * A
        v1 = 59.30 - 21.0 * (N - Z) / A - 0.024 * A
        v2 = 7.228e-3 - 1.48e-6 * A
        v3 = 1.994e-5 - 2.0e-8 * A
        w1 = 12.195 + 0.0167 * A
        d1 = 16.0 - 16.0 * (N - Z) / A
        avd = 0.5446 - 1.656e-4 * A
        awd = avd
        rc = 0.0
    elseif k0 == 2  # Protons
        ef = -8.4075 + 0.01378 * A
        v1 = 59.30 + 21.0 * (N - Z) / A - 0.024 * A
        v2 = 7.067e-3 + 4.23e-6 * A
        v3 = 1.729e-5 + 1.136e-8 * A
        w1 = 14.667 + 0.009629 * A
        avd = 0.5187 + 5.205e-4 * A
        awd = avd
        d1 = 16.0 + 16.0 * (N - Z) / A
        rc = 1.198 + 0.697 * A^(-2/3) + 12.994 * A^(-5/3)
    else
        error("Invalid particle type k0=$k0. Use 1 for neutron, 2 for proton.")
    end

    return (rv=rv, av=av, v1=v1, v2=v2, v3=v3, v4=v4, rw=rw, aw=aw, w1=w1, w2=w2,
            rvd=rvd, avd=avd, vd=vd, rwd=rwd, awd=awd, d1=d1, d2=d2, d3=d3,
            rvso=rvso, avso=avso, vso1=vso1, vso2=vso2, rwso=rwso, awso=awso,
            wso1=wso1, wso2=wso2, ef=ef, rc=rc)
end

"""
    energyform(k0, Z, A, Elab, params)

Calculate energy-dependent potential depths.

# Arguments
- `k0::Int`: Particle type (1 = neutron, 2 = proton)
- `Z::Float64`: Target charge number
- `A::Float64`: Target mass number
- `Elab::Float64`: Laboratory energy (MeV)
- `params`: Named tuple from globalomp

# Returns
Named tuple with (v, w, vd, wd, vso, wso) - potential depths in MeV.
"""
function energyform(k0::Int, Z::Float64, A::Float64, Elab::Float64, params)
    (; v1, v2, v3, v4, w1, w2, d1, d2, d3, vso1, vso2, wso1, wso2, ef, rc) = params

    f = Elab - ef

    # Coulomb correction
    if k0 == 1  # Neutron
        vcoul = 0.0
    else  # Proton
        Vc = 1.73 / rc * Z / A^(1/3)
        vcoul = Vc * v1 * (v2 - 2.0 * v3 * f + 3.0 * v4 * f^2)
    end

    # Potential depths
    v = v1 * (1.0 - v2 * f + v3 * f^2 - v4 * f^3) + vcoul
    w = w1 * f^2 / (f^2 + w2^2)
    vd = 0.0
    wd = d1 * f^2 * exp(-d2 * f) / (f^2 + d3^2)
    vso = vso1 * exp(-vso2 * f)
    wso = wso1 * f^2 / (f^2 + wso2^2)

    return (v=v, w=w, vd=vd, wd=wd, vso=vso, wso=wso)
end

"""
    kd02_potential(k0, Z, A, Elab; verbose=false)

Calculate Koning-Delaroche optical model potential parameters.

# Arguments
- `k0::Int`: Particle type (1 = neutron, 2 = proton)
- `Z::Float64`: Target charge number
- `A::Float64`: Target mass number
- `Elab::Float64`: Laboratory energy (MeV)
- `verbose::Bool`: If true, print potential parameters (default: false)

# Returns
`KD02Params` structure containing all optical model parameters.

# Example
```julia
# Proton on 208Pb at 30 MeV lab energy
params = kd02_potential(2, 82.0, 208.0, 30.0, verbose=true)
```

# Reference
A.J. Koning and J.P. Delaroche, Nucl. Phys. A713 (2003) 231-310
"""
function kd02_potential(k0::Int, Z::Float64, A::Float64, Elab::Float64; verbose::Bool=false)
    # Get global parameters
    gparams = globalomp(k0, Z, A)

    # Get energy-dependent depths
    edeps = energyform(k0, Z, A, Elab, gparams)

    # Create result structure
    result = KD02Params(
        edeps.v, gparams.rv, gparams.av,
        edeps.w, gparams.rw, gparams.aw,
        edeps.vd, gparams.rvd, gparams.avd,
        edeps.wd, gparams.rwd, gparams.awd,
        edeps.vso, gparams.rvso, gparams.avso,
        edeps.wso, gparams.rwso, gparams.awso,
        gparams.rc
    )

    if verbose
        parname = k0 == 1 ? "neutron" : "proton"
        println()
        println("Koning-Delaroche global optical model (June 2002)")
        println()
        println("          $parname on Z=$(Int(Z)) A=$(Int(A))")
        println()
        println("1. Optical model parameters, E=$(Elab) MeV")
        println()
        @printf("   V     rv    av     W     rw    aw\n")
        @printf("%7.3f%6.3f%6.3f%7.3f%6.3f%6.3f\n\n",
                result.v, result.rv, result.av, result.w, result.rw, result.aw)
        @printf("   Vd    rvd   avd    Wd    rwd   awd\n")
        @printf("%7.3f%6.3f%6.3f%7.3f%6.3f%6.3f\n\n",
                result.vd, result.rvd, result.avd, result.wd, result.rwd, result.awd)
        @printf("   Vso   rvso  avso   Wso   rwso  awso  rc\n")
        @printf("%7.3f%6.3f%6.3f%7.3f%6.3f%6.3f%6.3f\n\n",
                result.vso, result.rvso, result.avso, result.wso, result.rwso, result.awso, result.rc)
        println("----------------------------------------------")
    end

    return result
end

# Convenience method with integer arguments
kd02_potential(k0::Int, Z::Int, A::Int, Elab::Real; kwargs...) =
    kd02_potential(k0, Float64(Z), Float64(A), Float64(Elab); kwargs...)

end # module
