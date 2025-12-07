# SLAM.jl

**Scattering Lagrange Asymptotic Matching** - A Julia package for solving nuclear scattering problems using Lagrange mesh methods with Baye's exact differential matrices.

## Overview

SLAM.jl computes scattering cross-sections and S-matrix elements for nuclear reactions (e.g., proton-nucleus collisions) by solving the radial Schrödinger equation with Coulomb interactions. The package implements Baye's Lagrange-Legendre mesh method with exact differential matrices operating on expansion coefficients, providing high accuracy with minimal computational cost.

### Key Features

- **Lagrange-Legendre mesh** on (0, R) using Gauss-Legendre quadrature
- **Baye's exact D and T matrices** for analytically correct derivatives in coefficient space
- **x-regularized basis functions**: f_j(r) = α_j · r · L_j(r/R)
- **Coulomb wave functions** via COUL90 Fortran library (with pure Julia fallback)
- **Woods-Saxon optical model potentials** with volume, surface, and spin-orbit terms
- **Outgoing wave boundary conditions** using Coulomb-Hankel functions H⁺
- **Koning-Delaroche (KD02) global optical model** parameters
- **Numerov integration** for validation and comparison
- **Interactive GUI** for real-time calculations

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/your-username/SLAM.jl")
```

Or clone and develop locally:

```bash
git clone https://github.com/your-username/SLAM.jl
cd SLAM.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Dependencies

- Julia ≥ 1.6
- [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) - Gauss-Legendre quadrature
- [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl) - Gamma functions for phase shifts
- [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) - Matrix operations
- [Blink.jl](https://github.com/JuliaGizmos/Blink.jl) - Interactive GUI (optional)

## Quick Start

```julia
using SLAM

# Define optical potential for p + ⁴⁰Ca at 30 MeV
pot = OpticalPotential(
    V_v=50.0, r_v=1.25, a_v=0.65,    # Real volume
    W_v=10.0, r_wv=1.25, a_wv=0.65,  # Imaginary volume
    r_c=1.25,                         # Coulomb radius parameter
    Z_proj=1.0, A_proj=1.0,          # Proton
    Z_targ=20.0, A_targ=40.0         # ⁴⁰Ca
)

# Create scattering problem: E_cm = 30 MeV, l = 0
prob = ScatteringProblem(pot, 30.0, 0)

# Solve for S-matrix and cross section
result = solve_scattering(prob)

println("S-matrix:     ", result.S_matrix)
println("Phase shift:  ", result.phase_shift, " degrees")
println("Cross section:", result.cross_section, " mb")
```

## Usage Examples

### All Partial Waves

Sum over partial waves until convergence:

```julia
results = solve_all_partial_waves(pot, 30.0; l_max=20, N=60, R=20.0)
total_cs = sum(r.cross_section for r in results)
println("Total reaction cross section: $total_cs mb")
```

### Differential Cross Section

Compute angular distribution:

```julia
angles = 0:5:180  # degrees
dσ_dΩ = elastic_differential_cross_section(pot, 30.0, angles; l_max=25)
```

### Wave Function Extraction

Get the full scattering wave function:

```julia
wf = solve_wavefunction(prob)

# Access mesh points and wave function values
r_mesh = wf.r
ψ_scattered = wf.ψ_sc
ψ_total = wf.ψ_total

# Interpolate to arbitrary radius
ψ_at_5fm = evaluate_wavefunction(wf, 5.0)
```

### Numerov Comparison

Validate against standard Numerov outward integration:

```julia
num_result = solve_numerov(prob; h=0.05)

# Compare S-matrices
println("Lagrange S: ", result.S_matrix)
println("Numerov S:  ", num_result.S_matrix)
```

### Koning-Delaroche Global Potential

Use the KD02 energy-dependent global optical model:

```julia
include("src/kd02.jl")
using .KD02

# Get KD02 parameters for p + ¹²C at 30 MeV
params = kd02_potential(30.0, 12, 6, 2)  # E, A, Z, k0=2 for protons

# Convert to OpticalPotential
pot = OpticalPotential(
    V_v=params.v, r_v=params.rv, a_v=params.av,
    W_v=params.w, r_wv=params.rw, a_wv=params.aw,
    W_s=params.wd, r_ws=params.rwd, a_ws=params.awd,
    V_so=params.vso, r_so=params.rvso, a_so=params.avso,
    r_c=params.rc,
    Z_proj=1.0, A_proj=1.0,
    Z_targ=Float64(6), A_targ=Float64(12)
)
```

### Interactive GUI

Launch the graphical interface:

```bash
julia --project=. scripts/gui.jl
```

Or from Julia:

```julia
include("scripts/gui.jl")
```

## Mathematical Background

### Problem Formulation

SLAM.jl solves the radial Schrödinger equation in the center-of-mass frame:

$$\left[-\frac{\hbar^2}{2\mu}\frac{d^2}{dr^2} + \frac{\hbar^2 l(l+1)}{2\mu r^2} + V(r) - E\right]\psi(r) = 0$$

In dimensionless form (k² = 2μE/ℏ², U = 2μV/ℏ²):

$$\left[-\frac{d^2}{dr^2} + \frac{l(l+1)}{r^2} + U(r) - k^2\right]\psi(r) = 0$$

### Coulomb Scattering Decomposition

For charged-particle scattering, the wave function is decomposed:

$$\psi(r) = F_l(\eta, kr) + \psi_{\text{sc}}(r)$$

where F_l is the regular Coulomb function and ψ_sc satisfies the inhomogeneous equation:

$$\left[-\frac{d^2}{dr^2} + \frac{l(l+1)}{r^2} + U_{\text{full}} - k^2\right]\psi_{\text{sc}} = U_{\text{short}} \cdot F_l$$

### Lagrange-Legendre Mesh

The method uses N Gauss-Legendre quadrature points on (0,1), scaled to (0,R):
- Mesh points: r_j = R · x_j
- Quadrature weights: λ_j = w_j / 2
- x-regularized basis: f_j(r) = (r/r_j) · L_j(r/R) / √λ_j

The key relation connects expansion coefficients to function values:
$$c_j = \psi(r_j) \cdot \sqrt{R \cdot \lambda_j}$$

### Baye's Exact Matrices

The differential matrices act on coefficients rather than function values:

**D matrix** (d/dr scaled by R):
$$D_{ij} = (-1)^{i-j} \sqrt{\frac{x_i(1-x_j)}{x_j(1-x_i)}} \frac{1}{x_i - x_j} \quad (i \neq j)$$
$$D_{ii} = \frac{1}{2x_i(1-x_i)}$$

**T matrix** (-d²/dr² scaled by R²):
$$T_{ii} = \frac{N(N+1)x_i(1-x_i) - 3x_i + 1}{3x_i^2(1-x_i)^2}$$

### S-Matrix Extraction

From the boundary condition φ'(R) - γ_s·φ(R) = 0, where γ_s = k·H⁺'/H⁺:

$$S_l = 1 + 2ik \cdot f_l$$

with scattering amplitude f_l = ψ_sc(R) / [k · H⁺_l(kR)].

## Project Structure

```
SLAM.jl/
├── src/
│   ├── SLAM.jl           # Module root and exports
│   ├── mesh.jl           # Lagrange-Legendre mesh
│   ├── baye_matrices.jl  # Exact D and T matrices
│   ├── coulomb.jl        # Coulomb functions + COUL90 wrapper
│   ├── potentials.jl     # Woods-Saxon optical potentials
│   ├── scattering.jl     # Main scattering solver
│   ├── numerov.jl        # Numerov integration method
│   └── kd02.jl           # Koning-Delaroche global potential
├── test/
│   ├── runtests.jl       # Test suite
│   ├── compare_*.jl      # Method comparison scripts
│   └── debug_*.jl        # Diagnostic scripts
├── scripts/
│   └── gui.jl            # Interactive GUI
├── deps/
│   └── libcoul90.dylib   # Compiled COUL90 library
├── figures/              # Generated plots
└── Project.toml          # Package manifest
```

## API Reference

### Core Types

| Type | Description |
|------|-------------|
| `LagrangeMesh` | Mesh data: points, weights, normalization factors |
| `OpticalPotential` | Woods-Saxon parameters for nuclear + Coulomb potential |
| `ScatteringProblem` | Problem specification: potential, energy, angular momentum |
| `ScatteringResult` | Solution: S-matrix, phase shift, cross section |
| `WaveFunctionResult` | Wave function values on mesh + interpolation |
| `NumerovResult` | Numerov integration result for comparison |

### Core Functions

| Function | Description |
|----------|-------------|
| `init_legendre_mesh(N, R)` | Create N-point Lagrange-Legendre mesh on (0,R) |
| `baye_D_matrix(mesh)` | First derivative matrix |
| `baye_T_matrix(mesh)` | Kinetic energy matrix |
| `solve_scattering(prob)` | Solve single partial wave |
| `solve_all_partial_waves(pot, E)` | Sum all partial waves |
| `solve_wavefunction(prob)` | Get full wave function |
| `solve_numerov(prob)` | Numerov integration |
| `elastic_differential_cross_section(pot, E, θ)` | Angular distribution |

### Coulomb Functions

| Function | Description |
|----------|-------------|
| `coulomb_F(l, η, ρ)` | Regular Coulomb function F_l |
| `coulomb_G(l, η, ρ)` | Irregular Coulomb function G_l |
| `coulomb_FG(l, η, ρ)` | Both + derivatives |
| `coulomb_H_plus(l, η, ρ)` | Outgoing Hankel H⁺ = G + iF |
| `coulomb_eta(Z1, Z2, μ, E)` | Sommerfeld parameter |
| `coulomb_sigma(l, η)` | Coulomb phase shift |

## Testing

Run the test suite:

```julia
using Pkg
Pkg.test("SLAM")
```

Or directly:

```julia
include("test/runtests.jl")
```

Tests cover:
- Lagrange mesh initialization and properties
- Baye matrix structure
- Coulomb function Wronskian identity
- Optical potential evaluation
- Scattering solver convergence (|S| ≤ 1 for absorptive potentials)
- Legendre polynomial orthogonality

## Physical Constants

| Constant | Value | Description |
|----------|-------|-------------|
| ℏc | 197.3269718 MeV·fm | Reduced Planck constant × speed of light |
| e² | 1.43997 MeV·fm | Coulomb constant (includes 1/4πε₀) |
| AMU | 931.49432 MeV | Atomic mass unit |

## References

1. **D. Baye**, "The Lagrange-mesh method", *Physics Reports* 565 (2015) 1-107
   [doi:10.1016/j.physrep.2014.11.006](https://doi.org/10.1016/j.physrep.2014.11.006)

2. **A.R. Barnett**, "COULFG: Coulomb and Bessel functions and their derivatives", *Computer Physics Communications* 21 (1981) 297-314
   [doi:10.1016/0010-4655(81)90030-2](https://doi.org/10.1016/0010-4655(81)90030-2)

3. **A.J. Koning & J.P. Delaroche**, "Local and global nucleon optical models from 1 keV to 200 MeV", *Nuclear Physics A* 713 (2003) 231-310
   [doi:10.1016/S0375-9474(02)01321-0](https://doi.org/10.1016/S0375-9474(02)01321-0)

4. **NIST Digital Library of Mathematical Functions**, Chapter 33: Coulomb Functions
   [https://dlmf.nist.gov/33](https://dlmf.nist.gov/33)

## License

MIT License

## Author

Jin Lei (jinlei.phys@gmail.com)
